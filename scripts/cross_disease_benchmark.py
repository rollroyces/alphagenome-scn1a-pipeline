#!/usr/bin/env python3
"""
Cross-disease benchmark for AlphaGenome — 5 rare disease genes.

For each gene:
1. Extract SCN1A-style benchmark (216 pathogenic splicing + 375 benign intronic)
   - Adaptive: scales with available ClinVar pathogenic variants per gene
2. Score via AlphaGenome live API (SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS)
3. Compute AUROC, AUPRC, top-5% precision
4. Compare to SCN1A

Outputs:
- outputs/cross_disease_benchmark.csv  (one row per gene × scorer)
- outputs/cross_disease_results.md     (summary table + interpretation)
- figures/cross_disease_auprc.png      (bar chart)
"""

from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
)

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pysam


GENES = {
    "SCN1A":  {"chrom": "2",   "start": 165984640, "end": 166182806,
               "strand": "-", "transcript": "ENST00000674923.1", "disease": "Dravet syndrome"},
    "SCN2A":  {"chrom": "2",   "start": 165294829, "end": 165376368,
               "strand": "+", "transcript": "ENST00000283256.10", "disease": "Epileptic encephalopathy"},
    "MECP2":  {"chrom": "X",   "start": 154021572, "end": 154097731,
               "strand": "-", "transcript": "ENST00000303391.11", "disease": "Rett syndrome"},
    "CFTR":   {"chrom": "7",   "start": 117480025, "end": 117668665,
               "strand": "+", "transcript": "ENST00000003084.11", "disease": "Cystic fibrosis"},
    "DMD":    {"chrom": "X",   "start":  31119213, "end":  33339609,
               "strand": "-", "transcript": "ENST00000357033.9",  "disease": "Duchenne muscular dystrophy"},
}


def load_clinvar_for_gene(gene_name: str, info: dict, vcf_path: str, tabix_path: str,
                          n_positive: int = 200, n_negative: int = 350,
                          seed: int = 42) -> pd.DataFrame:
    """Extract pathogenic splicing + benign intronic variants for a gene from local ClinVar VCF."""
    tb = pysam.TabixFile(vcf_path, index=tabix_path)
    region = f"{info['chrom']}:{info['start']}-{info['end']}"
    rows = []
    try:
        for rec in tb.fetch(region):
            rows.append({
                "chrom": rec.contig, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0] if rec.alts else "",
                "info": rec.info,
            })
    except ValueError:
        # Try with chr prefix
        for rec in tb.fetch(f"chr{info['chrom']}:{info['start']}-{info['end']}"):
            rows.append({
                "chrom": rec.contig, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0] if rec.alts else "",
                "info": rec.info,
            })
    tb.close()

    # Parse info field
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df
    df["gene"] = df["info"].apply(lambda d: _get_info(d, "GENEINFO"))
    df["consequence"] = df["info"].apply(lambda d: _get_consequence(d))
    df["clnsig"] = df["info"].apply(lambda d: _get_info(d, "CLNSIG"))
    df["clnsig_category"] = df["clnsig"].apply(_classify_clnsig)
    df["gene"] = df["gene"].fillna("")
    df["gene_match"] = df["gene"].str.contains(gene_name, na=False)

    # Positive: pathogenic + splicing-related + gene match + SNV
    pos_mask = (
        df["gene_match"]
        & df["clnsig_category"].isin(["pathogenic"])
        & df["consequence"].str.contains("splice", na=False)
        & (df["ref"].str.len() == 1)
        & (df["alt"].str.len() == 1)
    )
    positives = df[pos_mask].copy()

    # Negative: benign + intronic + gene match + SNV
    neg_mask = (
        df["gene_match"]
        & df["clnsig_category"].isin(["benign"])
        & df["consequence"].str.contains("intron", na=False)
        & (df["ref"].str.len() == 1)
        & (df["alt"].str.len() == 1)
    )
    negatives = df[neg_mask].copy()

    # Subsample
    np.random.seed(seed)
    if len(positives) > n_positive:
        positives = positives.sample(n=n_positive, random_state=seed)
    if len(negatives) > n_negative:
        negatives = negatives.sample(n=n_negative, random_state=seed)

    out = pd.concat([positives, negatives], ignore_index=True)
    out["label"] = out.index.isin(positives.index).astype(int) if len(positives) > 0 else 0
    # Cleaner label
    out = out.reset_index(drop=True)
    out["label"] = "negative"
    out.loc[:len(positives)-1, "label"] = "positive"  # positives come first
    out["gene_name"] = gene_name
    print(f"  {gene_name}: {len(positives)} pathogenic splicing, {len(negatives)} benign intronic")
    return out


def _get_info(info: dict, key: str) -> str:
    val = info.get(key)
    if val is None:
        return ""
    if isinstance(val, str):
        return val
    if isinstance(val, (list, tuple)):
        return str(val[0]) if len(val) > 0 else ""
    return str(val)


def _get_consequence(info: dict) -> str:
    """Get MC consequence from ClinVar info dict."""
    mc = info.get("MC")
    if mc is None:
        return ""
    if isinstance(mc, (list, tuple)):
        return ",".join(str(x) for x in mc)
    return str(mc)


def _classify_clnsig(clnsig: str) -> str:
    if clnsig is None or clnsig == "":
        return "unknown"
    clnsig_lower = str(clnsig).lower()
    if "pathogenic" in clnsig_lower and "conflicting" not in clnsig_lower:
        if "likely" in clnsig_lower:
            return "likely_pathogenic"
        return "pathogenic"
    if "benign" in clnsig_lower and "conflicting" not in clnsig_lower:
        if "likely" in clnsig_lower:
            return "likely_benign"
        return "benign"
    if "uncertain" in clnsig_lower:
        return "uncertain"
    if "conflicting" in clnsig_lower:
        return "conflicting"
    return "other"


def score_variants(df: pd.DataFrame, gene_name: str, api_key: str) -> pd.DataFrame:
    """Score each variant with SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS."""
    dna_model = dna_client.create(api_key)
    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]
    scorer_names = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]

    rows = []
    t0 = time.time()
    n_success = 0
    n_fail = 0

    for i, row in df.iterrows():
        if i % 25 == 0 or i == len(df) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(df) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(df)}] {elapsed:.0f}s elapsed, "
                  f"success={n_success}, fail={n_fail}, ETA {remaining:.0f}s", flush=True)
        try:
            chrom = str(row["chrom"])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            variant = genome.Variant(
                chromosome=chrom,
                position=int(row["pos"]),
                reference_bases=str(row["ref"]),
                alternate_bases=str(row["alt"]),
            )
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
            scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
            result = {
                "gene": gene_name,
                "chrom": row["chrom"], "pos": row["pos"],
                "ref": row["ref"], "alt": row["alt"],
                "label": row["label"],
                "consequence": row.get("consequence", ""),
            }
            for ann, name in zip(scores, scorer_names):
                result[f"{name}_score"] = float(ann.X.sum())
            result["score_success"] = True
            n_success += 1
        except Exception as e:
            result = {
                "gene": gene_name,
                "chrom": row["chrom"], "pos": row["pos"],
                "ref": row["ref"], "alt": row["alt"],
                "label": row["label"],
                "score_success": False,
                "score_error": f"{type(e).__name__}: {str(e)[:80]}",
            }
            n_fail += 1
        rows.append(result)

    elapsed = time.time() - t0
    print(f"  Done. {n_success}/{len(df)} succeeded, {n_fail} failed, total {elapsed:.0f}s")
    return pd.DataFrame(rows)


def compute_metrics(df: pd.DataFrame, scorer_name: str = "SPLICE_SITES") -> dict:
    """Compute AUROC, AUPRC, top-5% precision for one scorer."""
    successful = df[df["score_success"] == True].copy()
    if len(successful) == 0 or successful["label"].nunique() < 2:
        return {"gene": df["gene"].iloc[0], "scorer": scorer_name,
                "n_total": len(successful), "n_positive": 0,
                "auroc": None, "auprc": None, "top_5_pct_precision": None}

    y_true = (successful["label"] == "positive").astype(int).values
    y_score = successful[f"{scorer_name}_score"].values

    auroc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))

    # Top-5% precision
    n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
    top_idx = np.argsort(y_score)[-n_top:]
    top_5_pct = float(y_true[top_idx].mean())

    return {
        "gene": df["gene"].iloc[0],
        "scorer": scorer_name,
        "n_total": len(successful),
        "n_positive": int(y_true.sum()),
        "n_negative": int((1 - y_true).sum()),
        "auroc": auroc,
        "auprc": auprc,
        "top_5_pct_precision": top_5_pct,
    }


def main() -> int:
    print("=" * 70)
    print("Cross-Disease AlphaGenome Benchmark — 5 rare disease genes")
    print("=" * 70)

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Use scripts/_run_with_key.sh")
        return 1

    vcf_path = "data/clinvar_grch38.vcf.gz"
    tabix_path = "data/clinvar_grch38.vcf.gz.tbi"

    if not os.path.exists(vcf_path):
        print(f"ERROR: ClinVar VCF not found at {vcf_path}")
        print(f"  Download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/")
        return 1

    all_results = []
    all_scores = []

    for gene_name, info in GENES.items():
        print(f"\n--- {gene_name} ({info['disease']}) ---")
        print(f"  Region: chr{info['chrom']}:{info['start']}-{info['end']}")

        # Extract ClinVar variants
        cv = load_clinvar_for_gene(gene_name, info, vcf_path, tabix_path)
        if len(cv) == 0 or cv["label"].nunique() < 2:
            print(f"  Skipping {gene_name}: insufficient benchmark data")
            all_results.append({"gene": gene_name, "scorer": "SPLICE_SITES",
                              "n_total": 0, "n_positive": 0,
                              "auroc": None, "auprc": None, "top_5_pct_precision": None,
                              "note": "insufficient data"})
            continue

        # Score
        scored = score_variants(cv, gene_name, api_key)
        scored.to_csv(f"outputs/cross_disease_{gene_name}_raw.csv", index=False)

        # Compute metrics
        for scorer_name in ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]:
            metrics = compute_metrics(scored, scorer_name)
            all_results.append(metrics)
            if metrics["auprc"] is not None:
                print(f"  {scorer_name}: AUPRC={metrics['auprc']:.4f}, "
                      f"AUROC={metrics['auroc']:.4f}, top-5%={metrics['top_5_pct_precision']:.3f}")

        all_scores.append(scored)

    # Save aggregate metrics
    os.makedirs("outputs", exist_ok=True)
    metrics_df = pd.DataFrame(all_results)
    metrics_df.to_csv("outputs/cross_disease_benchmark.csv", index=False)
    print(f"\nSaved → outputs/cross_disease_benchmark.csv")

    # Save all raw scores
    if all_scores:
        all_scores_df = pd.concat(all_scores, ignore_index=True)
        all_scores_df.to_csv("outputs/cross_disease_all_raw.csv", index=False)
        print(f"Saved → outputs/cross_disease_all_raw.csv ({len(all_scores_df)} rows)")

    # Summary table
    print("\n" + "=" * 70)
    print("CROSS-DISEASE SUMMARY (SPLICE_SITES)")
    print("=" * 70)
    splice_only = metrics_df[metrics_df["scorer"] == "SPLICE_SITES"].copy()
    print(splice_only[["gene", "n_total", "n_positive", "auroc", "auprc", "top_5_pct_precision"]].to_string(index=False))
    print("\n" + "=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
