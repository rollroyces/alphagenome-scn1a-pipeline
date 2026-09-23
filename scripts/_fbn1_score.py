#!/usr/bin/env python3
"""
Experiment 015 — Score FBN1 stratified variants with AlphaGenome live API.

Reads outputs/_fbn1_stratified.tsv (produced by _fbn1_extract.py),
builds two benchmark pools, and runs them sequentially:

  1. PRIMARY (UNFILTERED): positives = ALL pathogenic, negatives = ALL benign
                            (capped at 200 / 350)
  2. APPLES-TO-APPLES (FILTERED): positives = pathogenic + splice-related,
                            negatives = benign + intronic
                            (capped at 100 / 200)

FBN1 — chr15:48,408,312-48,645,721 (GENCODE v46, hg38, ~237 kb)
       ENSG00000166147, MANE Select ENST00000316623.10 (FBN1-201, minus strand)
       NCBI Gene ID: 2200
       Largest gene in the cross-disease benchmark (65 exons); first
       connective-tissue / fibroblast tissue representation.

Saves:
  outputs/cross_disease_fbn1_raw.csv
  outputs/cross_disease_fbn1_metrics.csv
  outputs/cross_disease_fbn1_filtered_raw.csv
  outputs/cross_disease_fbn1_filtered_metrics.csv

Run via:
    bash scripts/_run_with_key.sh scripts/_fbn1_score.py
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


SCORER_NAMES = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]

# Unfiltered caps (mirror KCNQ2 / COL4A5 / SCN1A pattern: 200/350)
N_POS_CAP_UNFILT = 200
N_NEG_CAP_UNFILT = 350

# Filtered caps: 323 pathogenic+splice and 861 benign+intronic available;
# cap to keep within the brief's 30-100 pos, 100-300 neg target
N_POS_CAP_FILT = 100
N_NEG_CAP_FILT = 200

SEED = 42


def make_unfiltered_benchmark(input_tsv: str) -> pd.DataFrame:
    df = pd.read_csv(input_tsv, sep="\t")
    pos = df[df["clnsig_category"] == "pathogenic"].copy()
    neg = df[df["clnsig_category"] == "benign"].copy()

    print(f"  pathogenic available: {len(pos)}, capped at {N_POS_CAP_UNFILT}")
    print(f"  benign available: {len(neg)}, capped at {N_NEG_CAP_UNFILT}")

    np.random.seed(SEED)
    if len(pos) > N_POS_CAP_UNFILT:
        pos = pos.sample(n=N_POS_CAP_UNFILT, random_state=SEED)
    if len(neg) > N_NEG_CAP_UNFILT:
        neg = neg.sample(n=N_NEG_CAP_UNFILT, random_state=SEED)

    pos["label"] = "positive"
    neg["label"] = "negative"
    bench = pd.concat([pos, neg], ignore_index=True)
    bench = bench.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    print(f"  Final unfiltered benchmark: {len(bench)} variants "
          f"({(bench['label']=='positive').sum()} pos, "
          f"{(bench['label']=='negative').sum()} neg)")
    return bench


def make_filtered_benchmark(input_tsv: str) -> pd.DataFrame:
    df = pd.read_csv(input_tsv, sep="\t")
    df["mc_l"] = df["mc"].fillna("").str.lower()
    pos = df[(df["clnsig_category"] == "pathogenic")
             & (df["mc_l"].str.contains("splice_donor_variant")
                | df["mc_l"].str.contains("splice_acceptor_variant")
                | df["mc_l"].str.contains("splice_region_variant"))].copy()
    neg = df[(df["clnsig_category"] == "benign")
             & (df["mc_l"].str.contains("intron_variant"))].copy()

    print(f"  pathogenic+splice available: {len(pos)}, capped at {N_POS_CAP_FILT}")
    print(f"  benign+intronic available: {len(neg)}, capped at {N_NEG_CAP_FILT}")

    np.random.seed(SEED)
    if len(pos) > N_POS_CAP_FILT:
        pos = pos.sample(n=N_POS_CAP_FILT, random_state=SEED)
    if len(neg) > N_NEG_CAP_FILT:
        neg = neg.sample(n=N_NEG_CAP_FILT, random_state=SEED)

    pos["label"] = "positive"
    neg["label"] = "negative"
    bench = pd.concat([pos, neg], ignore_index=True)
    bench = bench.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    print(f"  Final filtered benchmark: {len(bench)} variants "
          f"({(bench['label']=='positive').sum()} pos, "
          f"{(bench['label']=='negative').sum()} neg)")
    return bench


def score_one(dna_model, chrom: str, pos: int, ref: str, alt: str) -> dict:
    try:
        chrom_str = str(chrom)
        if not chrom_str.startswith("chr"):
            chrom_str = "chr" + chrom_str
        variant = genome.Variant(
            chromosome=chrom_str,
            position=int(pos),
            reference_bases=str(ref),
            alternate_bases=str(alt),
        )
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
        scorers = [
            variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in SCORER_NAMES
        ]
        scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
        )
        out = {"score_success": True}
        for ann, name in zip(scores, SCORER_NAMES):
            out[f"{name}_score"] = float(ann.X.sum())
        return out
    except Exception as e:
        return {
            "score_success": False,
            "score_error": f"{type(e).__name__}: {str(e)[:120]}",
        }


def bootstrap_auprc_ci(y_true, y_score, n_boot: int = 1000,
                        seed: int = 42) -> tuple[float, float]:
    """Bootstrap 95% CI on AUPRC (percentile method)."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    aucs = []
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(average_precision_score(y_true[idx], y_score[idx]))
    if not aucs:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def compute_metrics(df: pd.DataFrame, scorer_name: str, gene: str) -> dict:
    successful = df[df["score_success"] == True].copy()
    if len(successful) == 0 or successful["label"].nunique() < 2:
        return {
            "gene": gene, "scorer": scorer_name,
            "n_total": int(len(successful)), "n_positive": 0, "n_negative": 0,
            "auroc": None, "auprc": None,
            "auprc_ci_lo": None, "auprc_ci_hi": None,
            "top_5_pct_precision": None,
        }

    y_true = (successful["label"] == "positive").astype(int).values
    y_score = successful[f"{scorer_name}_score"].values

    auroc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))
    ci_lo, ci_hi = bootstrap_auprc_ci(y_true, y_score)

    n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
    top_idx = np.argsort(y_score)[-n_top:]
    top_5 = float(y_true[top_idx].mean())

    return {
        "gene": gene,
        "scorer": scorer_name,
        "n_total": int(len(successful)),
        "n_positive": int(y_true.sum()),
        "n_negative": int((1 - y_true).sum()),
        "auroc": auroc,
        "auprc": auprc,
        "auprc_ci_lo": ci_lo,
        "auprc_ci_hi": ci_hi,
        "top_5_pct_precision": top_5,
    }


def run_one(dna_model, bench: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    print(f"\n--- {label} run: scoring {len(bench)} variants ---")
    rows = []
    t0 = time.time()
    n_success = 0
    n_fail = 0
    n = len(bench)

    for i, row in bench.iterrows():
        if i % 25 == 0 or i == n - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (n - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{n}] {elapsed:.0f}s elapsed, "
                  f"success={n_success}, fail={n_fail}, ETA {remaining:.0f}s",
                  flush=True)
        result = score_one(dna_model, row["chrom"], row["pos"],
                           row["ref"], row["alt"])
        result.update({
            "chrom": row["chrom"], "pos": row["pos"],
            "ref": row["ref"], "alt": row["alt"],
            "rsid": row.get("rsid", ""),
            "label": row["label"],
            "clnsig_category": row.get("clnsig_category", ""),
            "clnsig_raw": row.get("clnsig", ""),
        })
        rows.append(result)
        if result["score_success"]:
            n_success += 1
        else:
            n_fail += 1

    elapsed = time.time() - t0
    print(f"\nDone ({label}). {n_success}/{n} succeeded, {n_fail} failed, total {elapsed:.0f}s")

    raw_df = pd.DataFrame(rows)

    metrics = []
    for scorer in SCORER_NAMES:
        m = compute_metrics(raw_df, scorer, gene="FBN1")
        metrics.append(m)
        if m["auprc"] is not None:
            print(f"  {scorer}: AUPRC={m['auprc']:.4f} "
                  f"[{m['auprc_ci_lo']:.3f}, {m['auprc_ci_hi']:.4f}], "
                  f"AUROC={m['auroc']:.4f}, top-5%={m['top_5_pct_precision']:.3f}, "
                  f"n={m['n_total']}")
    metrics_df = pd.DataFrame(metrics)
    return raw_df, metrics_df


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="outputs/_fbn1_stratified.tsv")
    parser.add_argument("--raw-out", default="outputs/cross_disease_fbn1_raw.csv")
    parser.add_argument("--metrics-out",
                        default="outputs/cross_disease_fbn1_metrics.csv")
    parser.add_argument("--filtered-raw-out",
                        default="outputs/cross_disease_fbn1_filtered_raw.csv")
    parser.add_argument("--filtered-metrics-out",
                        default="outputs/cross_disease_fbn1_filtered_metrics.csv")
    parser.add_argument("--skip-unfiltered", action="store_true")
    parser.add_argument("--skip-filtered", action="store_true")
    args = parser.parse_args()

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Use scripts/_run_with_key.sh")
        return 1

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    # Smoke test
    print("\nSmoke test...")
    try:
        v = genome.Variant(chromosome="chr15", position=48_500_000,
                            reference_bases="A", alternate_bases="G")
        interval = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
        _ = dna_model.score_variant(
            interval=interval, variant=v,
            variant_scorers=[
                variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"]
            ],
        )
        print("  smoke OK")
    except Exception as e:
        print(f"  smoke FAILED: {type(e).__name__}: {e}")
        return 2

    if not args.skip_unfiltered:
        print("\nBuilding FBN1 unfiltered benchmark...")
        unfiltered_bench = make_unfiltered_benchmark(args.input)
        raw_df, metrics_df = run_one(dna_model, unfiltered_bench, "unfiltered")
        raw_df.to_csv(args.raw_out, index=False)
        print(f"Saved → {args.raw_out}")
        metrics_df.to_csv(args.metrics_out, index=False)
        print(f"Saved → {args.metrics_out}")

    if not args.skip_filtered:
        print("\nBuilding FBN1 filtered (apples-to-apples) benchmark...")
        filtered_bench = make_filtered_benchmark(args.input)
        raw_df, metrics_df = run_one(dna_model, filtered_bench, "filtered")
        raw_df.to_csv(args.filtered_raw_out, index=False)
        print(f"Saved → {args.filtered_raw_out}")
        metrics_df.to_csv(args.filtered_metrics_out, index=False)
        print(f"Saved → {args.filtered_metrics_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())