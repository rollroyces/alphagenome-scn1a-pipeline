#!/usr/bin/env python3
"""
Benchmark AlphaGenome splicing predictions against ClinVar pathogenicity labels.

This is the head-line analysis. Given:
- outputs/clinvar_scn1a.tsv (labeled: pathogenic / benign / uncertain / etc.)
- outputs/atlas_scn1a_summary.csv (AlphaGenome splicing scores per variant)

It produces:
- AUPRC per splicing modality (the main number)
- Calibration curves (do scores correspond to probability of pathogenicity?)
- Top-K precision / recall curves
- A merged dataset saved to outputs/benchmark_scn1a.csv
- Figures in figures/

The headline number: how well does AlphaGenome's SPLICE_JUNCTIONS score
separate pathogenic from benign variants in SCN1A?

If AUPRC > 0.5 → AlphaGenome is informative for SCN1A splicing variants
If AUPRC > 0.7 → strong enough to publish as a benchmark paper
If AUPRC < 0.3 → model is not useful for this task; reframe the paper

Usage:
    # First run scripts/clinvar_scn1a_local.py and scripts/atlas_scn1a.py
    python scripts/benchmark_scn1a.py
"""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


# Splicing score columns we'll evaluate
SPLICING_SCORE_COLUMNS = [
    "SPLICE_SITES_score",
    "SPLICE_SITE_USAGE_score",
    "SPLICE_JUNCTIONS_score",
    "SPLICE_JUNCTIONS_ACTIVE_score",
]

SPLICING_QUANTILE_COLUMNS = [
    "SPLICE_SITES_quantile",
    "SPLICE_SITE_USAGE_quantile",
    "SPLICE_JUNCTIONS_quantile",
    "SPLICE_JUNCTIONS_ACTIVE_quantile",
]


@dataclass
class BenchmarkResult:
    """Per-scorer benchmark results."""

    scorer_name: str
    score_column: str
    n_variants_scored: int
    n_pathogenic: int
    n_benign: int
    n_uncertain: int
    auroc: float
    auprc: float
    top_1pct_precision: float
    top_5pct_precision: float


def load_clinvar(path: str) -> pd.DataFrame:
    """Load ClinVar SCN1A variants with categorical labels."""
    df = pd.read_csv(path, sep="\t")
    print(f"Loaded {len(df):,} ClinVar variants from {path}")
    return df


def load_atlas(path: str) -> pd.DataFrame:
    """Load AlphaGenome Atlas summary CSV."""
    df = pd.read_csv(path)
    print(f"Loaded {len(df):,} Atlas variants from {path}")
    return df


def normalize_variant_key(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Build a consistent (chrom, pos, ref, alt) key for joining.

    Args:
      source: 'clinvar' or 'atlas' — different column names.
    """
    df = df.copy()
    if source == "clinvar":
        # ClinVar uses numeric chrom (e.g., '2'), we keep as string
        df["chrom"] = df["chrom"].astype(str)
        df["pos"] = df["pos"].astype(int)
        df["ref"] = df["ref"].astype(str).str.upper()
        df["alt"] = df["alt"].astype(str).str.upper()
    elif source == "atlas":
        # Atlas uses 'chr2' format — normalize to '2'
        df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)
        df["pos"] = df["pos"].astype(int)
        df["ref"] = df["ref"].astype(str).str.upper()
        df["alt"] = df["alt"].astype(str).str.upper()
    return df


def join_clinvar_atlas(clinvar: pd.DataFrame, atlas: pd.DataFrame) -> pd.DataFrame:
    """Inner-join ClinVar with Atlas scores on (chrom, pos, ref, alt)."""
    joined = clinvar.merge(
        atlas,
        on=["chrom", "pos", "ref", "alt"],
        how="inner",
        suffixes=("_clinvar", "_atlas"),
    )
    print(f"Joined: {len(clinvar):,} ClinVar × {len(atlas):,} Atlas → {len(joined):,} matched variants")
    return joined


def compute_benchmark(joined: pd.DataFrame, score_col: str) -> BenchmarkResult | None:
    """Compute AUROC, AUPRC, top-K precision for one scorer."""
    # Filter to rows where this score is available
    scored = joined[joined[score_col].notna()].copy()
    if len(scored) == 0:
        return None

    # Build binary pathogenic label
    # Positive class: pathogenic / likely_pathogenic
    # Negative class: benign / likely_benign
    # We exclude uncertain_significance and conflicting from the AUROC computation
    binary = scored[scored["clnsig_category"].isin(["pathogenic", "benign"])].copy()
    binary["is_pathogenic"] = (binary["clnsig_category"] == "pathogenic").astype(int)

    n_scored = len(scored)
    n_pathogenic = int(binary["is_pathogenic"].sum())
    n_benign = int((1 - binary["is_pathogenic"]).sum())
    n_uncertain = int((scored["clnsig_category"] == "uncertain").sum())

    if n_pathogenic < 5 or n_benign < 5:
        print(f"  ⚠ {score_col}: too few labeled variants ({n_pathogenic} P / {n_benign} B), skipping")
        return None

    scores = binary[score_col].values
    labels = binary["is_pathogenic"].values

    # AlphaGenome scores are "disruption" — higher = more disruptive.
    # For pathogenicity, higher disruption → more likely pathogenic.
    # So no sign flip needed; treat raw score as pathogenicity score.

    try:
        auroc = roc_auc_score(labels, scores)
    except ValueError:
        auroc = float("nan")

    try:
        auprc = average_precision_score(labels, scores)
    except ValueError:
        auprc = float("nan")

    # Top-K precision: among variants with the top K% of scores, what fraction are pathogenic?
    sorted_scores = np.sort(scores)[::-1]
    n_top_1pct = max(1, int(0.01 * len(sorted_scores)))
    n_top_5pct = max(1, int(0.05 * len(sorted_scores)))

    threshold_1pct = sorted_scores[n_top_1pct - 1]
    threshold_5pct = sorted_scores[n_top_5pct - 1]

    top_1pct_labels = labels[scores >= threshold_1pct]
    top_5pct_labels = labels[scores >= threshold_5pct]

    top_1pct_precision = float(top_1pct_labels.mean()) if len(top_1pct_labels) > 0 else 0.0
    top_5pct_precision = float(top_5pct_labels.mean()) if len(top_5pct_labels) > 0 else 0.0

    return BenchmarkResult(
        scorer_name=score_col.replace("_score", "").replace("_ACTIVE", "_ACTIVE"),
        score_column=score_col,
        n_variants_scored=n_scored,
        n_pathogenic=n_pathogenic,
        n_benign=n_benign,
        n_uncertain=n_uncertain,
        auroc=float(auroc),
        auprc=float(auprc),
        top_1pct_precision=top_1pct_precision,
        top_5pct_precision=top_5pct_precision,
    )


def plot_pr_curves(joined: pd.DataFrame, results: list[BenchmarkResult], out_path: str):
    """Plot precision-recall curves for all splicing scorers."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for r in results:
        binary = joined[joined["clnsig_category"].isin(["pathogenic", "benign"])].copy()
        binary = binary[binary[r.score_column].notna()]
        binary["is_pathogenic"] = (binary["clnsig_category"] == "pathogenic").astype(int)

        precision, recall, _ = precision_recall_curve(
            binary["is_pathogenic"].values,
            binary[r.score_column].values,
        )
        ax.plot(recall, precision, label=f"{r.scorer_name} (AUPRC={r.auprc:.3f})", linewidth=2)

    baseline = joined[joined["clnsig_category"].isin(["pathogenic", "benign"])]["clnsig_category"].eq("pathogenic").mean()
    ax.axhline(baseline, color="gray", linestyle="--", alpha=0.5, label=f"baseline = {baseline:.3f}")

    ax.set_xlabel("Recall (pathogenic variants captured)")
    ax.set_ylabel("Precision (fraction of variants that are pathogenic)")
    ax.set_title("AlphaGenome splicing scores predict SCN1A pathogenicity\n(ClinVar pathogenic vs benign)")
    ax.legend(loc="lower left", fontsize=9)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.05])
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out_path}")


def plot_roc_curves(joined: pd.DataFrame, results: list[BenchmarkResult], out_path: str):
    """Plot ROC curves for all splicing scorers."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for r in results:
        binary = joined[joined["clnsig_category"].isin(["pathogenic", "benign"])].copy()
        binary = binary[binary[r.score_column].notna()]
        binary["is_pathogenic"] = (binary["clnsig_category"] == "pathogenic").astype(int)

        fpr, tpr, _ = roc_curve(
            binary["is_pathogenic"].values,
            binary[r.score_column].values,
        )
        ax.plot(fpr, tpr, label=f"{r.scorer_name} (AUROC={r.auroc:.3f})", linewidth=2)

    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", alpha=0.5, label="random")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate (sensitivity)")
    ax.set_title("AlphaGenome splicing scores vs SCN1A pathogenicity (ROC)")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.05])
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out_path}")


def plot_calibration(joined: pd.DataFrame, results: list[BenchmarkResult], out_path: str):
    """Plot calibration curves: predicted quantile vs observed pathogenicity rate.

    For each decile of predicted score, what's the actual fraction pathogenic?
    A well-calibrated model lies on the diagonal.
    """
    fig, axes = plt.subplots(1, len(results), figsize=(5 * len(results), 4), sharey=True)
    if len(results) == 1:
        axes = [axes]

    for ax, r in zip(axes, results):
        binary = joined[joined["clnsig_category"].isin(["pathogenic", "benign"])].copy()
        binary = binary[binary[r.score_column].notna()]

        # Use quantile column if available, else rank-transform scores
        qcol = r.score_column.replace("_score", "_quantile")
        if qcol in binary.columns and binary[qcol].notna().any():
            scores = binary[qcol].values
            xlabel = "Predicted quantile (calibrated)"
        else:
            scores = binary[r.score_column].values
            xlabel = "Raw score (rank-transformed)"

        labels = binary["clnsig_category"].eq("pathogenic").astype(int).values

        # Bin into deciles
        n_bins = 10
        sorted_idx = np.argsort(scores)
        bin_size = len(scores) // n_bins
        bin_pred = []
        bin_obs = []
        bin_count = []
        for i in range(n_bins):
            start = i * bin_size
            end = (i + 1) * bin_size if i < n_bins - 1 else len(scores)
            bin_idx = sorted_idx[start:end]
            bin_scores = scores[bin_idx]
            bin_labels = labels[bin_idx]
            bin_pred.append(bin_scores.mean())
            bin_obs.append(bin_labels.mean())
            bin_count.append(len(bin_idx))

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="perfect calibration")
        ax.scatter(bin_pred, bin_obs, s=[c * 5 for c in bin_count], alpha=0.7, color="C0")
        ax.plot(bin_pred, bin_obs, "C0-", alpha=0.5)
        ax.set_xlabel(xlabel)
        ax.set_title(r.scorer_name)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1.05])
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Observed pathogenicity rate")
    fig.suptitle("Calibration: AlphaGenome splicing scores vs SCN1A pathogenicity")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out_path}")


def plot_score_distribution(joined: pd.DataFrame, score_col: str, out_path: str):
    """Plot histogram of scores by clinical significance category."""
    fig, ax = plt.subplots(figsize=(8, 5))
    categories = ["pathogenic", "uncertain", "benign", "conflicting"]
    colors = {"pathogenic": "C3", "uncertain": "C1", "benign": "C2", "conflicting": "C0"}

    for cat in categories:
        subset = joined[joined["clnsig_category"] == cat]
        subset = subset[subset[score_col].notna()]
        if len(subset) > 0:
            ax.hist(subset[score_col], bins=50, alpha=0.5, label=f"{cat} (n={len(subset)})", color=colors[cat])

    ax.set_xlabel(f"{score_col}")
    ax.set_ylabel("Number of variants")
    ax.set_title(f"Score distribution: {score_col}")
    ax.legend()
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out_path}")


def main() -> int:
    clinvar_path = "outputs/clinvar_scn1a.tsv"
    atlas_path = "outputs/atlas_scn1a_summary.csv"
    if not os.path.exists(clinvar_path):
        print(f"ERROR: {clinvar_path} not found. Run scripts/clinvar_scn1a_local.py first.")
        return 1
    if not os.path.exists(atlas_path):
        print(f"ERROR: {atlas_path} not found. Run scripts/atlas_scn1a.py first.")
        return 1

    print("=" * 70)
    print("Loading data...")
    print("=" * 70)
    clinvar = load_clinvar(clinvar_path)
    atlas = load_atlas(atlas_path)

    clinvar = normalize_variant_key(clinvar, "clinvar")
    atlas = normalize_variant_key(atlas, "atlas")

    print("\n" + "=" * 70)
    print("Joining ClinVar with Atlas scores...")
    print("=" * 70)
    joined = join_clinvar_atlas(clinvar, atlas)
    print(f"After join: {len(joined):,} variants with both ClinVar label and Atlas score")
    print(f"Clinical significance of joined variants:")
    print(joined["clnsig_category"].value_counts().to_string())

    # Save merged dataset
    os.makedirs("outputs", exist_ok=True)
    merged_path = "outputs/benchmark_scn1a.csv"
    joined.to_csv(merged_path, index=False)
    print(f"Saved merged dataset → {merged_path}")

    print("\n" + "=" * 70)
    print("Running benchmark per splicing scorer...")
    print("=" * 70)

    results: list[BenchmarkResult] = []
    available_score_cols = [c for c in SPLICING_SCORE_COLUMNS if c in joined.columns]
    if not available_score_cols:
        print("ERROR: no splicing score columns found in Atlas summary")
        print(f"  Looked for: {SPLICING_SCORE_COLUMNS}")
        print(f"  Atlas columns: {list(atlas.columns)}")
        return 1
    print(f"Found score columns: {available_score_cols}")

    for col in available_score_cols:
        result = compute_benchmark(joined, col)
        if result is not None:
            results.append(result)
            print(f"\n  {result.scorer_name}:")
            print(f"    Scored: {result.n_variants_scored:,}")
            print(f"    Pathogenic: {result.n_pathogenic:,} | Benign: {result.n_benign:,} | Uncertain: {result.n_uncertain:,}")
            print(f"    AUROC: {result.auroc:.4f}")
            print(f"    AUPRC: {result.auprc:.4f}  ← headline number")
            print(f"    Top 1% precision: {result.top_1pct_precision:.3f}")
            print(f"    Top 5% precision: {result.top_5pct_precision:.3f}")

    # Save results table
    results_path = "outputs/benchmark_scn1a_results.csv"
    with open(results_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "scorer_name", "score_column", "n_variants_scored",
            "n_pathogenic", "n_benign", "n_uncertain",
            "auroc", "auprc", "top_1pct_precision", "top_5pct_precision",
        ])
        writer.writeheader()
        for r in results:
            writer.writerow(r.__dict__)
    print(f"\nSaved → {results_path}")

    # Generate figures
    print("\n" + "=" * 70)
    print("Generating figures...")
    print("=" * 70)
    os.makedirs("figures", exist_ok=True)

    if results:
        plot_pr_curves(joined, results, "figures/benchmark_pr_curves.png")
        plot_roc_curves(joined, results, "figures/benchmark_roc_curves.png")
        plot_calibration(joined, results, "figures/benchmark_calibration.png")

        # Score distribution for the best scorer (highest AUPRC)
        best = max(results, key=lambda r: r.auprc)
        plot_score_distribution(joined, best.score_column, f"figures/benchmark_score_distribution_{best.scorer_name}.png")

        print(f"\nBest scorer: {best.scorer_name} (AUPRC = {best.auprc:.4f})")
        print(f"  Top 1% of variants by this score: {best.top_1pct_precision * 100:.1f}% are pathogenic")
        print(f"  Top 5% of variants by this score: {best.top_5pct_precision * 100:.1f}% are pathogenic")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Outputs:")
    print(f"  {merged_path}  ← merged dataset")
    print(f"  {results_path}  ← benchmark metrics")
    print(f"  figures/benchmark_pr_curves.png")
    print(f"  figures/benchmark_roc_curves.png")
    print(f"  figures/benchmark_calibration.png")
    print(f"  figures/benchmark_score_distribution_*.png")
    print(f"\nNext:")
    print(f"  - If AUPRC > 0.5: AlphaGenome is informative, proceed to candidate identification")
    print(f"  - If AUPRC > 0.7: publishable as methods paper")
    print(f"  - If AUPRC < 0.3: model is not useful for splicing in SCN1A; reframe")
    return 0


if __name__ == "__main__":
    sys.exit(main())
