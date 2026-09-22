#!/usr/bin/env python3
"""
Run AlphaGenome SPLICE_SITES scoring across multiple rare disease gene benchmarks.

Computes AUROC, AUPRC, top-5% precision per gene.
Also performs bootstrap confidence intervals for genes with N≥20 positives.

Outputs:
- outputs/cross_disease_metrics.csv      (per-gene × scorer metrics)
- outputs/cross_disease_summary.md        (human-readable report)
- figures/cross_disease_auprc.png        (bar chart)
"""

from __future__ import annotations

import os
import sys
import glob
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
)

warnings.filterwarnings("ignore")


SCORERS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]


def compute_metrics(df: pd.DataFrame, scorer: str) -> dict:
    """Compute AUROC, AUPRC, top-5% precision for one scorer, with bootstrap CI on AUPRC."""
    successful = df[df["score_success"] == True].copy()
    if len(successful) == 0 or successful["label"].nunique() < 2:
        return {"gene": df["gene"].iloc[0] if len(df) > 0 else "?",
                "scorer": scorer, "n_total": len(successful), "n_positive": 0,
                "auroc": None, "auprc": None, "auprc_ci_lo": None, "auprc_ci_hi": None,
                "top_5_pct_precision": None,
                "note": "insufficient data"}

    y_true = (successful["label"] == "positive").astype(int).values
    y_score = successful[f"{scorer}_score"].values

    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())

    auroc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))

    # Top-5% precision
    n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
    top_idx = np.argsort(y_score)[-n_top:]
    top_5_pct = float(y_true[top_idx].mean())

    # Bootstrap CI on AUPRC (only if n_pos >= 10 — otherwise too noisy)
    rng = np.random.default_rng(42)
    n_boot = 200
    boot_auprcs = []
    if n_pos >= 10:
        pos_idx = np.where(y_true == 1)[0]
        neg_idx = np.where(y_true == 0)[0]
        for _ in range(n_boot):
            # Stratified bootstrap
            bs_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
            bs_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
            bs_idx = np.concatenate([bs_pos, bs_neg])
            bs_y = y_true[bs_idx]
            bs_score = y_score[bs_idx]
            if bs_y.sum() == 0 or bs_y.sum() == len(bs_y):
                continue
            boot_auprcs.append(average_precision_score(bs_y, bs_score))
    if len(boot_auprcs) >= 50:
        ci_lo = float(np.percentile(boot_auprcs, 2.5))
        ci_hi = float(np.percentile(boot_auprcs, 97.5))
    else:
        ci_lo = ci_hi = None

    return {
        "gene": df["gene"].iloc[0],
        "scorer": scorer,
        "n_total": len(successful),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "auroc": auroc,
        "auprc": auprc,
        "auprc_ci_lo": ci_lo,
        "auprc_ci_hi": ci_hi,
        "top_5_pct_precision": top_5_pct,
    }


def load_all_results(pattern: str = "outputs/cross_disease_*_raw.csv") -> pd.DataFrame:
    """Load and concatenate all per-gene scoring results."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files match {pattern}")
    print(f"Loading {len(files)} per-gene scoring files:")
    dfs = []
    for f in files:
        df = pd.read_csv(f)
        print(f"  {f}: {len(df)} rows")
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def plot_results(metrics_df: pd.DataFrame, output_path: str):
    """Generate bar chart of AUPRC per gene."""
    splice = metrics_df[(metrics_df["scorer"] == "SPLICE_SITES")
                        & metrics_df["auprc"].notna()].copy()
    splice = splice.sort_values("auprc", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar with confidence intervals where available
    x_pos = np.arange(len(splice))
    bars = ax.bar(x_pos, splice["auprc"],
                  yerr=[splice["auprc"] - splice["auprc_ci_lo"].fillna(splice["auprc"]),
                        splice["auprc_ci_hi"].fillna(splice["auprc"]) - splice["auprc"]],
                  capsize=4, alpha=0.7, color="steelblue", edgecolor="navy")

    # Color top performer differently
    if len(splice) > 0:
        bars[0].set_color("darkorange")
        bars[0].set_edgecolor("darkred")

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(splice["gene"], splice["n_positive"])],
                       fontsize=11)
    ax.set_ylabel("AUPRC (AlphaGenome SPLICE_SITES)", fontsize=11)
    ax.set_title("Cross-disease AlphaGenome benchmark\n"
                 "Pathogenic splicing vs. benign intronic variant discrimination",
                 fontsize=13)
    ax.set_ylim(0, 1.05)
    ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5, label="AUPRC=0.5 (random)")
    ax.grid(axis="y", alpha=0.3)

    # Annotate values on bars
    for i, (g, a, n) in enumerate(zip(splice["gene"], splice["auprc"], splice["n_positive"])):
        ax.text(i, a + 0.02, f"{a:.3f}", ha="center", fontsize=10, fontweight="bold")
        ax.text(i, -0.05, f"n_pos={n}", ha="center", fontsize=8, color="gray")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved → {output_path}")


def write_summary(metrics_df: pd.DataFrame, output_path: str):
    """Write human-readable summary."""
    with open(output_path, "w") as f:
        f.write("# Cross-Disease AlphaGenome Benchmark — Results\n\n")
        f.write("SPLICE_SITES scores from AlphaGenome for pathogenic splicing variants vs. benign intronic variants\n")
        f.write("across 5 rare disease genes. Higher AUPRC = better discrimination.\n\n")

        # Per-scorer summary
        for scorer in SCORERS:
            f.write(f"## {scorer}\n\n")
            sub = metrics_df[metrics_df["scorer"] == scorer].copy()
            sub = sub.sort_values("auprc", ascending=False, na_position="last")
            f.write("| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |\n")
            f.write("|------|---------|-------|-------|-------|-------|--------|-------------|\n")
            for _, row in sub.iterrows():
                if row["auprc"] is None or pd.isna(row["auprc"]):
                    f.write(f"| {row['gene']} | {row['n_total']} | 0 | 0 | — | — | — | — |\n")
                else:
                    ci = ""
                    if row["auprc_ci_lo"] is not None and not pd.isna(row["auprc_ci_lo"]):
                        ci = f"[{row['auprc_ci_lo']:.3f}, {row['auprc_ci_hi']:.3f}]"
                    f.write(f"| {row['gene']} | {row['n_total']} | {row['n_positive']} | "
                           f"{row['n_negative']} | {row['auroc']:.4f} | {row['auprc']:.4f} | "
                           f"{ci} | {row['top_5_pct_precision']:.3f} |\n")
            f.write("\n")

        # Interpretation
        f.write("## Interpretation\n\n")
        splice = metrics_df[(metrics_df["scorer"] == "SPLICE_SITES")
                           & metrics_df["auprc"].notna()].sort_values("auprc", ascending=False)
        if len(splice) > 0:
            f.write(f"- **Best performer:** {splice.iloc[0]['gene']} (AUPRC={splice.iloc[0]['auprc']:.4f})\n")
            f.write(f"- **Worst performer:** {splice.iloc[-1]['gene']} (AUPRC={splice.iloc[-1]['auprc']:.4f})\n")
            f.write(f"- **Mean AUPRC across genes:** {splice['auprc'].mean():.4f}\n")
            f.write(f"- **Std AUPRC across genes:** {splice['auprc'].std():.4f}\n\n")

            # Notes
            f.write("### Notes\n\n")
            for _, row in splice.iterrows():
                notes = []
                if row["n_positive"] < 20:
                    notes.append("small positive set, confidence intervals are wide")
                if row["auprc_ci_lo"] is not None and row["auprc_ci_hi"] is not None:
                    width = row["auprc_ci_hi"] - row["auprc_ci_lo"]
                    if width < 0.05:
                        notes.append("tight CI")
                    elif width > 0.15:
                        notes.append("wide CI")
                if notes:
                    f.write(f"- **{row['gene']}** (n_pos={row['n_positive']}): {', '.join(notes)}\n")

    print(f"Saved → {output_path}")


def main() -> int:
    print("=" * 70)
    print("Cross-Disease Benchmark Aggregator")
    print("=" * 70)

    # Load all per-gene results
    df = load_all_results()

    # Compute metrics
    all_metrics = []
    for gene in df["gene"].unique():
        gene_df = df[df["gene"] == gene]
        for scorer in SCORERS:
            metrics = compute_metrics(gene_df, scorer)
            all_metrics.append(metrics)

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv("outputs/cross_disease_metrics.csv", index=False)
    print(f"\nSaved → outputs/cross_disease_metrics.csv")

    # Plot
    os.makedirs("figures", exist_ok=True)
    plot_results(metrics_df, "figures/cross_disease_auprc.png")

    # Summary
    write_summary(metrics_df, "outputs/cross_disease_summary.md")

    # Print to console
    print("\n" + "=" * 70)
    print("CROSS-DISEASE SUMMARY (SPLICE_SITES)")
    print("=" * 70)
    splice = metrics_df[(metrics_df["scorer"] == "SPLICE_SITES")
                        & metrics_df["auprc"].notna()].sort_values("auprc", ascending=False)
    print(f"{'Gene':<8} {'n_pos':<8} {'AUROC':<8} {'AUPRC':<8} {'95% CI':<20} {'Top-5%':<8}")
    print("-" * 60)
    for _, row in splice.iterrows():
        ci = ""
        if row["auprc_ci_lo"] is not None and not pd.isna(row["auprc_ci_lo"]):
            ci = f"[{row['auprc_ci_lo']:.3f}, {row['auprc_ci_hi']:.3f}]"
        else:
            ci = "n/a (small N)"
        print(f"{row['gene']:<8} {row['n_positive']:<8} {row['auroc']:.4f}   "
              f"{row['auprc']:.4f}   {ci:<20} {row['top_5_pct_precision']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
