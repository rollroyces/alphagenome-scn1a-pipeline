#!/usr/bin/env python3
"""
Benchmark AlphaGenome SCN1A splicing predictions against ClinVar — using live API.

The Atlas SDK has a thread-pool bug that fails at ~47% of large queries. As a
workaround, we score each variant individually with the live prediction API
(score_variant), which doesn't have that bug.

For 216 pathogenic splicing + ~500 benign intronic controls in SCN1A, that's
~700 API calls. At ~1/sec, that's ~12 minutes — within the free tier limit.

Output: outputs/benchmark_scn1a.csv (merged dataset with AlphaGenome scores)
        figures/benchmark_*.png (PR, ROC, calibration, score distribution)

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/benchmark_scn1a_live_api.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# SCN1A locus
SCN1A_CHROM = "chr2"  # dna_client uses 'chr' prefix
SCN1A_START = 165_984_640
SCN1A_END = 166_182_806
CONTEXT_BP = 100_000  # 100 Kb context each side — but we use 16 Kb windows for cost


def load_clinvar() -> pd.DataFrame:
    path = "outputs/clinvar_scn1a.tsv"
    df = pd.read_csv(path, sep="\t")
    # Build positive controls: pathogenic + splicing-related
    pos = df[(df["clnsig_category"] == "pathogenic") & (df["is_splicing_related"] == "yes")].copy()
    print(f"Loaded {len(df):,} total ClinVar SCN1A variants")
    print(f"  Positive controls (pathogenic + splicing): {len(pos):,}")
    # Build negative controls: benign intronic
    neg = df[(df["clnsig_category"] == "benign") & (df["molecular_consequence"] == "intron_variant")].copy()
    print(f"  Negative controls (benign intronic): {len(neg):,}")
    # Sample negatives to keep API calls manageable — 3x positives
    if len(neg) > 3 * len(pos):
        neg = neg.sample(n=3 * len(pos), random_state=42).copy()
        print(f"  Down-sampled negatives to {len(neg):,}")
    return pd.concat([pos, neg], ignore_index=True)


def score_one_variant(dna_model: dna_client.DnaClient, variant_row: pd.Series) -> dict:
    """Score one variant with all 3 splicing scorers via the live API."""
    try:
        # Build variant — note: dna_client uses 'chr' prefix in Interval, but
        # variant.chromosome is the bare chromosome (no 'chr')? Actually it does
        # accept 'chr2' in the constructor.
        variant = genome.Variant(
            chromosome=("chr" + str(variant_row["chrom"])) if not str(variant_row["chrom"]).startswith("chr") else str(variant_row["chrom"]),
            position=int(variant_row["pos"]),
            reference_bases=str(variant_row["ref"]),
            alternate_bases=str(variant_row["alt"]),
        )
        # 16 Kb window centered on variant (cheapest option)
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)

        # Use SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
        scorers_to_use = [
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
        ]

        scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers_to_use,
        )
        # scores is a list of AnnData (one per scorer)
        result = {"success": True}
        for ann, scorer_name in zip(scores, ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]):
            score_sum = float(ann.X.sum())
            result[f"{scorer_name}_score"] = score_sum
        return result
    except Exception as e:
        return {"success": False, "error": f"{type(e).__name__}: {e}"}


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    print("Loading ClinVar controls...")
    controls = load_clinvar()
    print(f"\nTotal variants to score: {len(controls):,}")
    print(f"Estimated API calls: {len(controls)} (each variant is one call)")

    print("\nCreating dna_model client...")
    dna_model = dna_client.create(api_key)

    output_rows = []
    t0 = time.time()

    for i, row in controls.iterrows():
        if i % 25 == 0 or i == len(controls) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(controls) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(controls)}] {elapsed:.0f}s elapsed, ETA {remaining:.0f}s")

        score_result = score_one_variant(dna_model, row)

        merged_row = dict(row)
        merged_row.update(score_result)
        output_rows.append(merged_row)

    elapsed = time.time() - t0
    print(f"\nFinished scoring in {elapsed:.0f}s")

    # Save raw results
    raw_path = "outputs/benchmark_scn1a_live_api_raw.csv"
    out_df = pd.DataFrame(output_rows)
    out_df.to_csv(raw_path, index=False)
    print(f"Saved raw → {raw_path}")

    # Compute benchmark
    successful = out_df[out_df["success"] == True].copy()
    failed = out_df[out_df["success"] != True]
    print(f"\nSuccessful: {len(successful):,} / {len(out_df):,}")
    if len(failed) > 0:
        print(f"Failed: {len(failed):,}")
        print("First few failures:")
        print(failed[["chrom", "pos", "ref", "alt", "error"]].head(5).to_string())

    if len(successful) == 0:
        print("No successful scores — cannot benchmark")
        return 2

    # Add binary pathogenic label
    successful["is_pathogenic"] = (successful["clnsig_category"] == "pathogenic").astype(int)

    print("\n=== Benchmark: SCN1A AlphaGenome splicing scores vs ClinVar pathogenicity ===")
    from sklearn.metrics import (
        average_precision_score,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    score_cols = ["SPLICE_SITES_score", "SPLICE_SITE_USAGE_score", "SPLICE_JUNCTIONS_score"]

    results = []
    for col in score_cols:
        if col not in successful.columns:
            continue
        valid = successful[successful[col].notna()]
        if len(valid) < 5:
            continue
        y_true = valid["is_pathogenic"].values
        y_score = valid[col].values
        try:
            auroc = roc_auc_score(y_true, y_score)
            auprc = average_precision_score(y_true, y_score)
            n_pos = int(y_true.sum())
            n_neg = int((1 - y_true).sum())
            # Top-K precision
            sorted_scores = np.sort(y_score)[::-1]
            n_top5 = max(1, int(0.05 * len(sorted_scores)))
            threshold = sorted_scores[n_top5 - 1]
            top5_labels = y_true[y_score >= threshold]
            top5_prec = float(top5_labels.mean()) if len(top5_labels) > 0 else 0.0
            results.append({"scorer": col, "n": len(valid), "n_pos": n_pos, "n_neg": n_neg, "auroc": auroc, "auprc": auprc, "top5_precision": top5_prec})
            print(f"\n  {col}:")
            print(f"    N = {len(valid)} (pos={n_pos}, neg={n_neg})")
            print(f"    AUROC = {auroc:.4f}")
            print(f"    AUPRC = {auprc:.4f}  ← headline number")
            print(f"    Top 5% precision = {top5_prec:.3f}")
        except ValueError as e:
            print(f"  {col}: cannot compute — {e}")

    # Save results table
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv("outputs/benchmark_scn1a_results.csv", index=False)
        print(f"\nSaved → outputs/benchmark_scn1a_results.csv")

        # Generate PR curve figure
        os.makedirs("figures", exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 6))
        for r in results:
            valid = successful[successful[r["scorer"]].notna()]
            y_true = valid["is_pathogenic"].values
            y_score = valid[r["scorer"]].values
            prec, rec, _ = precision_recall_curve(y_true, y_score)
            ax.plot(rec, prec, label=f'{r["scorer"].replace("_score", "")} (AUPRC={r["auprc"]:.3f})', linewidth=2)
        baseline = successful["is_pathogenic"].mean()
        ax.axhline(baseline, color="gray", linestyle="--", alpha=0.5, label=f"baseline = {baseline:.3f}")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title(f"AlphaGenome vs SCN1A pathogenicity\n({int(baseline*100)}% baseline prevalence, {len(successful)} ClinVar variants)")
        ax.legend(loc="lower left")
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1.05])
        plt.tight_layout()
        plt.savefig("figures/benchmark_pr_curves_live_api.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved → figures/benchmark_pr_curves_live_api.png")

        # Score distribution plot
        fig, ax = plt.subplots(figsize=(10, 5))
        best_col = max(results, key=lambda r: r["auprc"])["scorer"]
        for cat, color in [("pathogenic", "C3"), ("uncertain", "C1"), ("benign", "C2")]:
            subset = successful[successful["clnsig_category"] == cat]
            subset = subset[subset[best_col].notna()]
            if len(subset) > 0:
                ax.hist(subset[best_col], bins=30, alpha=0.5, label=f"{cat} (n={len(subset)})", color=color)
        ax.set_xlabel(f"{best_col}")
        ax.set_ylabel("Number of variants")
        ax.set_title(f"Score distribution: {best_col}")
        ax.legend()
        plt.tight_layout()
        plt.savefig("figures/benchmark_score_dist_live_api.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved → figures/benchmark_score_dist_live_api.png")

        # Headline
        best = max(results, key=lambda r: r["auprc"])
        print(f"\n★ BEST SCORER: {best['scorer']}")
        print(f"  AUPRC = {best['auprc']:.4f}")
        print(f"  Top 5% precision = {best['top5_precision']:.3f}")
        if best["auprc"] > 0.7:
            print("  → Strong enough for methods paper submission")
        elif best["auprc"] > 0.5:
            print("  → Informative; publishable as benchmark")
        elif best["auprc"] > 0.3:
            print("  → Marginal; consider reframing")
        else:
            print("  → Below baseline; AlphaGenome not useful here")

    return 0


if __name__ == "__main__":
    sys.exit(main())
