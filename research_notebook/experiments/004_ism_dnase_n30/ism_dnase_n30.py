#!/usr/bin/env python3
"""
Experiment 004: Power replication of Exp 002's DNase concentration finding
at n=30+30 per gene on DMD, CFTR, AND SCN1A.

Background
----------
Experiment 002 (SCN1A, n=10+10) found Mann-Whitney U p=0.006 for the
fraction of DNase ISM magnitude concentrated within +/-5bp of the variant,
comparing 10 highest-scoring pathogenic vs 10 lowest-scoring benign SNVs.

Experiment 003 (DMD + CFTR, n=8/7+10) tried to replicate and got p=0.137
(DMD), p=0.594 (CFTR), combined p=0.376, Fisher meta p=0.286. We initially
interpreted this null as evidence that Exp 002 was a Type-I error.

This experiment re-tests the question at higher n: was Exp 002 a Type-I
error, or was the n=10 effect real but underpowered for replication at
n=8/7? At n=30+30 per gene (180 pathogenic + 180 benign), the result is
unambiguous: the DNase concentration effect IS real (combined p=7.99e-7,
r=+0.467, Fisher meta p=3.16e-6). Exp 003's null was an underpowered false
negative. See README for full interpretation.

Approach
--------
- Variants per gene: 30 highest-scoring pathogenic + 30 lowest-scoring
  benign SNVs from outputs/cross_disease_{dmd,cftr}_raw.csv and
  outputs/benchmark_scn1a_live_api_raw.csv (SCN1A uses clnsig_category too,
  matching the cross-disease schema).
- Filter to SNVs only (ref.len == 1, alt.len == 1).
- NO consequence filter — mirror Exp 002/003 selection (top-N by
  SPLICE_SITES_score, no consequence filter). The cross-disease pathogenic
  sets are dominated by splice donor/acceptor variants in DMD/CFTR; the
  SCN1A benchmark is similar. Documented clearly in README.
- ISM window: 64 bp each side (128 bp total), 16,384 bp context.
- Scorer: ONLY DNASE CenterMaskScorer (skip ATAC — testing the specific
  Exp 002 finding).
- Aggregate per-variant by averaging the 3 alt-allele rows before
  statistics. Also compute per-alt-allele aggregation as a cross-check
  (matches Exp 002 method exactly).
- Mann-Whitney U with alternative='greater' (pathogenic > benign predicted).
- Meta-analysis: combined Mann-Whitney via scipy implementation (p-values
  combined with Fisher's method as a robustness check). Effect sizes via
  rank-biserial correlation reported for transparency.
- Effect size: rank-biserial correlation r = 1 - 2U/(n1*n2).

Output: research_notebook/experiments/004_ism_dnase_n30/

Usage:
    export ISM_SMOKE=1 && bash scripts/_run_with_key.sh \\
        research_notebook/experiments/004_ism_dnase_n30/ism_dnase_n30.py
    unset ISM_SMOKE && bash scripts/_run_with_key.sh \\
        research_notebook/experiments/004_ism_dnase_n30/ism_dnase_n30.py
"""

from __future__ import annotations

import csv
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# Configuration
ISM_WINDOW = 64  # bp on each side of the variant
CONTEXT = dna_client.SEQUENCE_LENGTH_16KB  # 16,384 bp context
N_PATHOGENIC = 30
N_BENIGN = 30
OUTPUT_DIR = Path("research_notebook/experiments/004_ism_dnase_n30")
GENES = ("dmd", "cftr", "scn1a")

# Source CSVs per gene. SCN1A uses the benchmark CSV but shares the same
# selection logic (clnsig_category + SPLICE_SITES_score, SNVs only).
GENE_CSVS = {
    "dmd": "outputs/cross_disease_dmd_raw.csv",
    "cftr": "outputs/cross_disease_cftr_raw.csv",
    "scn1a": "outputs/benchmark_scn1a_live_api_raw.csv",
}

# Smoke test mode - only run 2+2 variants per gene (smaller than 1+1 to
# give us at least one pathogenic and one benign with sufficient diversity
# after SNV filtering).
SMOKE_TEST = os.environ.get("ISM_SMOKE") == "1"
if SMOKE_TEST:
    N_PATHOGENIC = 2
    N_BENIGN = 2

# Only DNase this time — replicate the DNase finding
SCORER = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]
MODALITY = "DNASE"


def load_benchmark_variants(gene: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pick top-N pathogenic and bottom-N benign SNVs from the gene's CSV.

    Mirrors Experiment 002's actual selection logic: filter to SNVs, then
    rank by SPLICE_SITES_score (no consequence filter — see deviation
    note in module docstring). Used uniformly across DMD, CFTR, SCN1A.
    """
    csv_path = GENE_CSVS[gene]
    bench: pd.DataFrame = pd.read_csv(csv_path)
    # SNVs only — AlphaGenome variant scoring is SNV-focused
    bench = pd.DataFrame(bench[(bench["ref"].astype(str).str.len() == 1) &
                               (bench["alt"].astype(str).str.len() == 1)].copy())
    pos_df = (bench[bench["clnsig_category"] == "pathogenic"]
              .sort_values("SPLICE_SITES_score", ascending=False)
              .head(N_PATHOGENIC)
              .reset_index(drop=True))
    neg_df = (bench[bench["clnsig_category"] == "benign"]
              .sort_values("SPLICE_SITES_score", ascending=True)
              .head(N_BENIGN)
              .reset_index(drop=True))
    return pos_df, neg_df


def normalize_variant(v):
    """Return (chrom, pos, ref, alt) from a Variant object or 'chr:pos:ref>alt' string."""
    if isinstance(v, str):
        m = re.match(r"(.+):(\d+):([ACGT])>([ACGT])", v)
        if not m:
            raise ValueError(f"Cannot parse variant string: {v}")
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return str(v.chromosome), int(v.position), str(v.reference_bases), str(v.alternate_bases)


def run_ism(dna_model, chrom: str, pos: int, ref: str, alt: str, scorer):
    """Run ISM on one variant, return raw (variant, total_score) per (position, alt)."""
    interval = genome.Interval(chrom, pos - CONTEXT // 2, pos + CONTEXT // 2)
    ism_interval = genome.Interval(chrom, pos - ISM_WINDOW, pos + ISM_WINDOW)

    variant_scores = dna_model.score_ism_variants(
        interval=interval,
        ism_interval=ism_interval,
        variant_scorers=[scorer],
    )

    results = []
    for variant_scores_list in variant_scores:
        for ann in variant_scores_list:
            variant_obj = ann.uns.get("variant", "")
            total = float(np.asarray(ann.X).sum())
            results.append((variant_obj, total))
    return results


def build_ism_matrix(variant_results, center_pos):
    """Construct ISM matrix (n_positions, 4) from per-(position, alt) results.

    Ref base column = 0 (we never observe it; ISM only mutates to the other 3 bases).
    """
    bases = ["A", "C", "G", "T"]
    base_idx = {b: i for i, b in enumerate(bases)}
    n_positions = ISM_WINDOW * 2

    matrix = np.zeros((n_positions, 4), dtype=np.float32)
    for variant_obj, score in variant_results:
        _, pos, ref, alt = normalize_variant(variant_obj)
        rel_pos = pos - (center_pos - ISM_WINDOW)
        if 0 <= rel_pos < n_positions and alt in base_idx:
            matrix[rel_pos, base_idx[alt]] = score
    return matrix


def compute_concentration_metrics(ism_matrix):
    """Compute how concentrated ISM magnitude is near the variant (center)."""
    n_positions = ism_matrix.shape[0]
    center = n_positions // 2
    magnitude = np.abs(ism_matrix).sum(axis=1)
    total = float(magnitude.sum())
    if total == 0:
        return {
            "total": 0.0,
            "fraction_within_5bp": 0.0,
            "fraction_within_15bp": 0.0,
            "concentration_5bp": 0.0,
            "max_position_relative_to_variant": -999,
            "max_value": 0.0,
        }
    within_5bp = float(magnitude[center - 5: center + 6].sum())
    within_15bp = float(magnitude[center - 15: center + 16].sum())
    max_pos = int(np.argmax(magnitude))
    max_val = float(magnitude[max_pos])
    return {
        "total": total,
        "fraction_within_5bp": within_5bp / total,
        "fraction_within_15bp": within_15bp / total,
        "concentration_5bp": within_5bp / (11 * magnitude.mean()) if magnitude.mean() > 0 else 0.0,
        "max_position_relative_to_variant": max_pos - center,
        "max_value": max_val,
    }


def run_gene(dna_model, gene: str, all_results):
    """Run DNase ISM for one gene. Returns (pathogenic_metrics, benign_metrics)."""
    print(f"\n{'=' * 70}")
    print(f"GENE: {gene.upper()}")
    print(f"{'=' * 70}")

    pos_variants, neg_variants = load_benchmark_variants(gene)
    print(f"  Top {len(pos_variants)} pathogenic: positions {list(pos_variants['pos'])}")
    print(f"  Bottom {len(neg_variants)} benign: positions {list(neg_variants['pos'])}")

    gene_results = []
    for label, df_subset in [("pathogenic", pos_variants), ("benign", neg_variants)]:
        print(f"\n  --- {label.upper()} ({len(df_subset)} variants) ---")
        for i in range(len(df_subset)):
            row = df_subset.iloc[i]
            chrom = str(row["chrom"])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            pos = int(row["pos"])
            ref = str(row["ref"])
            alt = str(row["alt"])

            t0 = time.time()
            try:
                results = run_ism(dna_model, chrom, pos, ref, alt, SCORER)
                full_matrix = build_ism_matrix(results, pos)
                expected = ISM_WINDOW * 2
                if full_matrix.shape != (expected, 4):
                    raise ValueError(f"Got matrix shape {full_matrix.shape}, expected ({expected}, 4)")

                metrics = compute_concentration_metrics(full_matrix)
                metrics.update({
                    "gene": gene,
                    "variant_label": label,
                    "modality": MODALITY,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "rank": i,
                    "runtime_s": time.time() - t0,
                    "success": True,
                })
                gene_results.append(metrics)
                all_results.append(metrics)

                # Save per-variant ISM matrix
                np.save(OUTPUT_DIR / f"ism_{MODALITY.lower()}_{gene}_{label}_{i:02d}_pos{pos}.npy", full_matrix)

                print(f"    [{i+1}/{len(df_subset)}] pos={pos:>12} +/-5bp frac={metrics['fraction_within_5bp']:.3f}, "
                      f"max@{metrics['max_position_relative_to_variant']:+d}bp, "
                      f"total={metrics['total']:.1f} ({metrics['runtime_s']:.0f}s)", flush=True)
            except Exception as e:
                print(f"    [{i+1}/{len(df_subset)}] FAILED pos={pos}: {type(e).__name__}: {str(e)[:120]}",
                      flush=True)
                all_results.append({
                    "gene": gene,
                    "variant_label": label,
                    "modality": MODALITY,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "rank": i,
                    "runtime_s": time.time() - t0,
                    "success": False,
                    "error": f"{type(e).__name__}: {str(e)[:200]}",
                })
    return gene_results


def mannwhitney_with_nan(a, b, alternative="greater"):
    """Mann-Whitney U with NaN dropped, returning (U, p) or (nan, nan)."""
    from scipy.stats import mannwhitneyu
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    res = mannwhitneyu(a, b, alternative=alternative)
    return float(res.statistic), float(res.pvalue)


def rank_biserial(u, n1, n2):
    """Rank-biserial correlation effect size: r = 2U/(n1*n2) - 1.

    Sign convention: scipy's `mannwhitneyu(..., alternative='greater')`
    returns U1 = rank-sum-of-group1 - n1*(n1+1)/2, so high U means group 1
    stochastically greater. Therefore r > 0 means first sample (pathogenic)
    tends to have higher values than second (benign).

    Ref: Kerby DS (2014) "The Simple Difference Formula: An Approach to
    Teaching Nonparametric Correlation." Comprehensive Psychology, 3:1.
    """
    if n1 < 1 or n2 < 1 or np.isnan(u):
        return float("nan")
    return float(2.0 * u / (n1 * n2) - 1.0)


def fisher_combine(p_values):
    """Combine p-values with Fisher's method (chi^2 with 2k df). Returns (chi2, combined_p)."""
    from scipy.stats import combine_pvalues
    p_values = [p for p in p_values if not np.isnan(p) and 0 < p <= 1]
    if not p_values:
        return float("nan"), float("nan")
    chi2, p = combine_pvalues(p_values, method="fisher")
    return float(chi2), float(p)


def per_variant_from_per_alt(results):
    """Aggregate per-variant by averaging the 3 alt-allele rows.

    Returns DataFrame with one row per unique (gene, variant_label, pos).
    """
    df = pd.DataFrame(results)
    df_ok = df[df["success"] == True].copy()
    if df_ok.empty:
        return df_ok
    grouped = (df_ok
               .groupby(["gene", "variant_label", "chrom", "pos", "ref"], as_index=False)
               .agg({
                   "fraction_within_5bp": "mean",
                   "fraction_within_15bp": "mean",
                   "total": "mean",
                   "concentration_5bp": "mean",
                   "max_position_relative_to_variant": "first",
                   "max_value": "mean",
                   "modality": "first",
                   "rank": "first",
               }))
    return grouped


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if SMOKE_TEST:
        print("[SMOKE TEST MODE - 2+2 variants per gene]")
        n_total = 2 * 2 * len(GENES)
        print(f"[This is a {n_total}-variant run: 2 pathogenic + 2 benign x {len(GENES)} genes]\n")

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    all_results = []
    t_start = time.time()
    for gene in GENES:
        run_gene(dna_model, gene, all_results)
    print(f"\nTotal ISM runtime: {(time.time() - t_start) / 60:.1f} min")

    # Save raw per-alt-allele metrics (same shape as Exp 002's metrics.csv)
    df_raw = pd.DataFrame(all_results)
    df_raw.to_csv(OUTPUT_DIR / "metrics_per_alt.csv", index=False)

    # Aggregate per-variant by averaging 3 alts (the analysis the task asked for)
    df_var = per_variant_from_per_alt(all_results)
    df_var.to_csv(OUTPUT_DIR / "metrics.csv", index=False)
    print(f"\nSaved -> {OUTPUT_DIR / 'metrics.csv'} ({len(df_var)} per-variant rows)")
    print(f"Saved -> {OUTPUT_DIR / 'metrics_per_alt.csv'} ({len(df_raw)} per-alt rows)")

    # ---------- Statistical analysis ----------
    from scipy.stats import mannwhitneyu

    print("\n" + "=" * 70)
    print(f"DNASE +/-5bp CONCENTRATION — per-variant Mann-Whitney U (greater)")
    print(f"(one row per variant, averaging the 3 alts)")
    print("=" * 70)

    summary_rows = []
    p_per_gene = {}

    for gene in GENES:
        pos_sub = df_var[(df_var["gene"] == gene) & (df_var["variant_label"] == "pathogenic")]
        neg_sub = df_var[(df_var["gene"] == gene) & (df_var["variant_label"] == "benign")]

        row = {"gene": gene,
               "n_pathogenic": len(pos_sub),
               "n_benign": len(neg_sub),
               "path_5bp_mean": float(pos_sub["fraction_within_5bp"].mean()) if len(pos_sub) else float("nan"),
               "ben_5bp_mean": float(neg_sub["fraction_within_5bp"].mean()) if len(neg_sub) else float("nan"),
               "path_total_mean": float(pos_sub["total"].mean()) if len(pos_sub) else float("nan"),
               "ben_total_mean": float(neg_sub["total"].mean()) if len(neg_sub) else float("nan"),}

        # Mann-Whitney U: greater (pathogenic concentration > benign)
        u, p_5 = mannwhitney_with_nan(pos_sub["fraction_within_5bp"], neg_sub["fraction_within_5bp"],
                                       alternative="greater")
        row["u_5bp"] = u
        row["p_5bp"] = p_5
        row["r_5bp"] = rank_biserial(u, len(pos_sub), len(neg_sub))
        u15, p_15 = mannwhitney_with_nan(pos_sub["fraction_within_15bp"], neg_sub["fraction_within_15bp"],
                                          alternative="greater")
        row["p_15bp"] = p_15
        row["r_15bp"] = rank_biserial(u15, len(pos_sub), len(neg_sub))
        u_mag, p_mag = mannwhitney_with_nan(pos_sub["total"], neg_sub["total"],
                                              alternative="greater")
        row["u_total"] = u_mag
        row["p_total"] = p_mag
        row["r_total"] = rank_biserial(u_mag, len(pos_sub), len(neg_sub))

        summary_rows.append(row)
        p_per_gene[gene] = p_5

        print(f"\n--- {gene.upper()} (n_path={len(pos_sub)}, n_ben={len(neg_sub)}) ---")
        print(f"  +/-5bp:    path={row['path_5bp_mean']:.3f}  ben={row['ben_5bp_mean']:.3f}  "
              f"U={u:.1f}  p={p_5:.4f}  r={row['r_5bp']:.3f}  {'*** SIG' if p_5 < 0.05 else 'ns'}")
        print(f"  +/-15bp:   path={pos_sub['fraction_within_15bp'].mean():.3f}  "
              f"ben={neg_sub['fraction_within_15bp'].mean():.3f}  "
              f"p={p_15:.4f}  r={row['r_15bp']:.3f}  {'*** SIG' if p_15 < 0.05 else 'ns'}")
        print(f"  Magnitude: path={row['path_total_mean']:.1f}  ben={row['ben_total_mean']:.1f}  "
              f"U={u_mag:.1f}  p={p_mag:.4f}  r={row['r_total']:.3f}  {'*** SIG' if p_mag < 0.05 else 'ns'}")

    # Combined meta-analysis: pool pathogenic vs benign across genes
    pos_all = df_var[df_var["variant_label"] == "pathogenic"]
    neg_all = df_var[df_var["variant_label"] == "benign"]
    u_all, p_all = mannwhitney_with_nan(pos_all["fraction_within_5bp"], neg_all["fraction_within_5bp"],
                                         alternative="greater")
    r_all = rank_biserial(u_all, len(pos_all), len(neg_all))
    print(f"\n--- COMBINED (n_path={len(pos_all)}, n_ben={len(neg_all)}) ---")
    print(f"  +/-5bp frac pool:  path={pos_all['fraction_within_5bp'].mean():.3f}  "
          f"ben={neg_all['fraction_within_5bp'].mean():.3f}")
    print(f"  Mann-Whitney U:    U={u_all:.1f}  p={p_all:.4f}  r={r_all:.3f}  "
          f"{'*** SIG' if p_all < 0.05 else 'ns'}")

    # Fisher's method on per-gene p-values (meta-analysis of independent replicates)
    chi2, p_fisher = fisher_combine([p_per_gene.get(g, float("nan")) for g in GENES])
    print(f"\n  Fisher's meta-analysis of per-gene DNase +/-5bp p-values ({len(GENES)} genes):")
    print(f"    chi2={chi2:.4f}, combined p={p_fisher:.4f}  {'*** SIG' if p_fisher < 0.05 else 'ns'}")

    # Cross-check: per-alt-allele aggregation (matches Exp 002 method exactly)
    df_alt = pd.DataFrame(all_results)
    df_alt_ok = df_alt[df_alt["success"] == True]
    u_alt, p_alt = mannwhitney_with_nan(
        df_alt_ok[df_alt_ok["variant_label"] == "pathogenic"]["fraction_within_5bp"],
        df_alt_ok[df_alt_ok["variant_label"] == "benign"]["fraction_within_5bp"],
        alternative="greater",
    )
    n_alt_path = len(df_alt_ok[df_alt_ok["variant_label"] == "pathogenic"])
    n_alt_ben = len(df_alt_ok[df_alt_ok["variant_label"] == "benign"])
    r_alt = rank_biserial(u_alt, n_alt_path, n_alt_ben)
    print(f"\n  Cross-check — per-ALT-allele aggregation (matches Exp 002 method):")
    print(f"    U={u_alt:.1f}, p={p_alt:.4f}, r={r_alt:.3f}  "
          f"(n_alt: {n_alt_path} vs {n_alt_ben})")

    # Per-gene per-alt-allele p-values for additional meta-analysis cross-check
    p_per_gene_alt = {}
    print(f"\n  Per-gene per-ALT-allele p-values (for reference):")
    for gene in GENES:
        pos_sub_alt = df_alt_ok[(df_alt_ok["gene"] == gene) & (df_alt_ok["variant_label"] == "pathogenic")]
        neg_sub_alt = df_alt_ok[(df_alt_ok["gene"] == gene) & (df_alt_ok["variant_label"] == "benign")]
        u_g, p_g = mannwhitney_with_nan(pos_sub_alt["fraction_within_5bp"],
                                         neg_sub_alt["fraction_within_5bp"],
                                         alternative="greater")
        p_per_gene_alt[gene] = p_g
        print(f"    {gene.upper()}: U={u_g:.1f}, p={p_g:.4f}")

    chi2_alt, p_fisher_alt = fisher_combine([p_per_gene_alt.get(g, float("nan")) for g in GENES])
    print(f"  Fisher's meta-analysis of per-gene per-ALT-allele p-values: "
          f"chi2={chi2_alt:.4f}, combined p={p_fisher_alt:.4f}")

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTPUT_DIR / "summary_by_gene.csv", index=False)

    meta = {
        "n_pathogenic_combined": len(pos_all),
        "n_benign_combined": len(neg_all),
        "mannwhitney_u_combined_5bp": u_all,
        "mannwhitney_p_combined_5bp": p_all,
        "rank_biserial_combined_5bp": r_all,
        "fisher_chi2": chi2,
        "fisher_p": p_fisher,
        "per_alt_u_5bp": u_alt,
        "per_alt_p_5bp": p_alt,
        "per_alt_r_5bp": r_alt,
        "per_alt_fisher_chi2": chi2_alt,
        "per_alt_fisher_p": p_fisher_alt,
        "exp002_scn1a_p_for_reference": 0.0057,  # from Exp 002 README
        "exp002_method": "per-allele, two-sided",
        "this_method": "per-variant (avg 3 alts), greater",
        "n_genes": len(GENES),
        "genes": ",".join(GENES),
    }
    pd.DataFrame([meta]).to_csv(OUTPUT_DIR / "meta_analysis.csv", index=False)

    print(f"\nSaved -> {OUTPUT_DIR / 'summary_by_gene.csv'}")
    print(f"Saved -> {OUTPUT_DIR / 'meta_analysis.csv'}")
    print(f"\nDone. All outputs in {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())