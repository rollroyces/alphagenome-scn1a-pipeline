#!/usr/bin/env python3
"""
Experiment 006: Tissue-specific vs averaged DNASE scoring for SCN1A (brain)
and DMD (muscle).

Background
----------
Experiment 004 (n=30+30 per gene, ISM concentration ±5bp) showed:
- DMD: p=7.97e-6, r=+0.695 (large effect, muscle gene)
- SCN1A: p=0.007, r=+0.421 (medium effect, brain gene)
- CFTR: p=0.318 (no effect, epithelial gene)

Experiment 005 (SPLICE_JUNCTIONS) showed tissue-filtered tracks did NOT
significantly outperform averaged tracks across 5 genes (Wilcoxon on
deltas, p>0.05).

Hypothesis
----------
Tissue-filtered DNASE tracks (brain for SCN1A, muscle for DMD) will give
a STRONGER pathogenic-vs-benign separation than averaged (all 305) DNASE
tracks, because AlphaGenome's tissue-specific track signal should be more
discriminative for genes whose disease phenotype is tissue-restricted.

Method
------
- Variants: SAME selection as Exp 004 — top-30 pathogenic + bottom-30
  benign SNVs per gene, ranked by SPLICE_SITES_score.
  - SCN1A: outputs/benchmark_scn1a_live_api_raw.csv
  - DMD: outputs/cross_disease_dmd_raw.csv
  - Excluded CFTR (task brief doesn't include it; we already have a
    null result for CFTR from Exp 004).
- Per-variant aggregation: average across 3 alt alleles, matching
  Exp 004's per-variant aggregation → final n=20+28 SCN1A, 24+29 DMD.

*** Important deviation from the task brief, documented honestly:
The task brief said "Run DNASE CenterMaskScorer ISM with 64bp window each
side (same as Exp 004)". However, `dna_client.score_ism_variants` does
NOT accept ontology_terms at the API level (only the `atlas.query_variants`
endpoint supports tissue filtering). Furthermore, the ISM AnnData in
Exp 004 is pre-aggregated across all 305 DNASE tracks BEFORE being
returned (the saved matrices are (128 positions, 4 bases) — track-level
resolution is lost).

So tissue filtering at the track level is impossible with ISM. We
therefore use `atlas.query_variants` (same endpoint as Exp 005) for DNASE
scoring, which returns AnnData with shape (n_variants, n_tracks=305) and
accepts `ontology_terms`. Per-variant score = max-abs across the
n_tracks dimension, as specified.

This is a deliberate deviation documented up-front so the comparison is
interpretable. The comparison still answers the substantive question:
"does tissue-filtered DNASE outperform averaged DNASE on pathogenic-vs-
benign separation?" — using identical variants and an identical
per-variant aggregation.
- Conditions:
  1. averaged: no `ontology_terms` → 305 tracks
  2. tissue-filtered:
     - SCN1A: brain-related DNASE tracks (biosample_name matches
       brain|neuron|cortex|hippocamp|cerebell|neural|glia|spinal|
       forebrain|midbrain, case-insensitive)
     - DMD: muscle-related DNASE tracks (muscle|myocyte|skeletal|
       cardiac|myoblast|smooth muscle, case-insensitive)
- Per-variant score: max(|X|) across tracks. (Secondary: mean(|X|).)
- Statistical test: Mann-Whitney U, alternative='greater' (pathogenic >
  benign predicted), rank-biserial r effect size.
- Output: per-variant scores per gene × condition, summary table,
  comparison plot, README with honest interpretation.

Usage:
    bash scripts/_run_with_key.sh \\
        research_notebook/experiments/006_tissue_specific_dnase/tissue_dnase_experiment.py
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from alphagenome.data import genome, ontology
from alphagenome.atlas import atlas

OUTPUT_DIR = Path("research_notebook/experiments/006_tissue_specific_dnase")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENES = ("scn1a", "dmd")

GENE_CSVS = {
    "scn1a": "outputs/benchmark_scn1a_live_api_raw.csv",
    "dmd": "outputs/cross_disease_dmd_raw.csv",
}

# Tissue keywords (case-insensitive biosample_name regex)
BRAIN_KEYWORDS = (
    "brain", "neuron", "cortex", "hippocamp", "cerebell", "neural",
    "glia", "spinal", "forebrain", "midbrain", "dorsolateral",
)
MUSCLE_KEYWORDS = (
    "muscle", "myocyte", "skeletal", "cardiac", "myoblast",
)
# Exclude 'smooth muscle of the brain vasculature' — that's a brain
# pericyte proxy, not skeletal/cardiac muscle. We'll filter it out.
MUSCLE_EXCLUDE = ("brain vasculature",)

# Per-gene tissue specification
GENE_TISSUE = {
    "scn1a": ("brain", BRAIN_KEYWORDS, ()),
    "dmd":    ("muscle", MUSCLE_KEYWORDS, MUSCLE_EXCLUDE),
}


def load_benchmark_variants(gene: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mirror Exp 004 selection: top-30 pathogenic + bottom-30 benign
    SNVs ranked by SPLICE_SITES_score."""
    bench = pd.read_csv(GENE_CSVS[gene])
    bench = bench[
        (bench["ref"].astype(str).str.len() == 1) &
        (bench["alt"].astype(str).str.len() == 1)
    ].copy()

    pos_df = (bench[bench["clnsig_category"] == "pathogenic"]
              .sort_values("SPLICE_SITES_score", ascending=False)
              .head(30)
              .reset_index(drop=True))
    neg_df = (bench[bench["clnsig_category"] == "benign"]
              .sort_values("SPLICE_SITES_score", ascending=True)
              .head(30)
              .reset_index(drop=True))
    return pos_df, neg_df


def make_variant(row) -> genome.Variant:
    chrom = str(row["chrom"]).strip()
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return genome.Variant(
        chromosome=chrom,
        position=int(row["pos"]),
        reference_bases=str(row["ref"]).strip(),
        alternate_bases=str(row["alt"]).strip(),
    )


def get_tissue_ontology_terms(atlas_client, tissue: str, keywords, exclude):
    """Return unique ontology terms matching tissue keywords in DNASE."""
    sm = atlas_client.scorer_metadata()
    dn = sm["DNASE"].track_metadata
    pattern_include = "|".join(keywords)
    mask = dn["biosample_name"].str.contains(pattern_include, case=False, na=False)
    if exclude:
        pattern_exclude = "|".join(exclude)
        mask &= ~dn["biosample_name"].str.contains(pattern_exclude, case=False, na=False)
    tracks = dn[mask]
    print(f"  {tissue.upper()} tracks matched: {len(tracks)} of {len(dn)}")
    print(f"  biosample_names: {sorted(tracks['biosample_name'].unique().tolist())}")

    terms = []
    for _, row in tracks.iterrows():
        curie = str(row["ontology_curie"])
        try:
            terms.append(ontology.from_curie(curie))
        except Exception:
            pass
    return list(set(terms)), tracks


def run_dnase(atlas_client, variants, ontology_terms=None, label=""):
    """Query DNASE for variants. Returns per-variant AnnData."""
    kwargs = dict(
        variants=variants,
        requested_scorers=["DNASE"],
        progress_bar=False,
        max_workers=4,
    )
    if ontology_terms is not None:
        kwargs["ontology_terms"] = ontology_terms

    t0 = time.time()
    result = atlas_client.query_variants(**kwargs)
    elapsed = time.time() - t0

    ad = result["DNASE"]
    print(f"  ({label}) shape={ad.shape}, elapsed={elapsed:.1f}s")
    return ad, elapsed


def anndata_to_per_variant_maxabs(ad, var_df) -> pd.DataFrame:
    """For each variant row in ad, compute max(|X[i, :]|) across tracks.

    Return DataFrame with columns: chrom, pos, ref, alt, max_abs, mean_abs,
    n_tracks. Merged by variant key (chr:pos:ref>alt) back to var_df.
    """
    X = ad.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    obs = ad.obs.copy()
    obs["max_abs"] = np.abs(X).max(axis=1)
    obs["mean_abs"] = np.abs(X).mean(axis=1)
    obs["n_tracks"] = X.shape[1]
    obs["variant_key"] = obs["variant"].astype(str)

    # Merge back to original variant info (variant_key == chr:pos:ref>alt)
    var_df = var_df.copy()

    def to_key(r):
        chrom = str(r["chrom"]).strip()
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        return f"{chrom}:{int(r['pos'])}:{r['ref']}>{r['alt']}"
    var_df["variant_key"] = var_df.apply(to_key, axis=1)

    merged = var_df.merge(
        obs[["variant_key", "max_abs", "mean_abs", "n_tracks"]],
        on="variant_key", how="left"
    )
    # Drop duplicate columns from merge
    if "variant_key" in merged.columns:
        merged = merged.drop(columns=["variant_key"])
    return merged


def mannwhitney_r(a, b, alternative="greater"):
    """Mann-Whitney U + rank-biserial r. NaN-safe.

    Sign convention: scipy's `mannwhitneyu(..., alternative='greater')`
    returns U1 = rank-sum-of-group1 - n1*(n1+1)/2, so high U means group 1
    stochastically greater. Therefore r = 2U/(n1*n2) - 1 (Kerby 2014):
    r > 0 means first sample (pathogenic) tends to have higher values
    than second (benign).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan"), float("nan"), len(a), len(b)
    res = mannwhitneyu(a, b, alternative=alternative)
    u = float(res.statistic)
    p = float(res.pvalue)
    n1, n2 = len(a), len(b)
    r = 2.0 * u / (n1 * n2) - 1.0
    return u, p, r, n1, n2


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    print("Creating atlas client...")
    atlas_client = atlas.create(api_key)

    # Pre-compute tissue ontology terms for each gene
    tissue_terms_by_gene = {}
    print("\n=== Tissue track discovery ===")
    for gene in GENES:
        tissue_label, keywords, exclude = GENE_TISSUE[gene]
        print(f"\n{gene.upper()} ({tissue_label}):")
        terms, tracks = get_tissue_ontology_terms(atlas_client, tissue_label,
                                                   keywords, exclude)
        tissue_terms_by_gene[gene] = terms

    all_per_variant = []  # for plotting
    summary_rows = []

    for gene in GENES:
        print(f"\n{'=' * 70}")
        print(f"GENE: {gene.upper()}")
        print(f"{'=' * 70}")
        pos_df, neg_df = load_benchmark_variants(gene)
        print(f"  Top-30 pathogenic: {len(pos_df)} rows, "
              f"{pos_df.groupby(['pos','ref']).ngroups} unique (pos,ref)")
        print(f"  Bottom-30 benign: {len(neg_df)} rows, "
              f"{neg_df.groupby(['pos','ref']).ngroups} unique (pos,ref)")

        pos_variants = [make_variant(r) for _, r in pos_df.iterrows()]
        neg_variants = [make_variant(r) for _, r in neg_df.iterrows()]

        for cond, ontology_terms, label in [
            ("avg", None, "averaged (305 tracks)"),
            ("tissue", tissue_terms_by_gene[gene], "tissue-filtered"),
        ]:
            print(f"\n  --- Condition: {label} ---")
            try:
                pos_ad, pos_elapsed = run_dnase(
                    atlas_client, pos_variants,
                    ontology_terms=ontology_terms,
                    label=f"pathogenic/{cond}",
                )
                neg_ad, neg_elapsed = run_dnase(
                    atlas_client, neg_variants,
                    ontology_terms=ontology_terms,
                    label=f"benign/{cond}",
                )
            except Exception as e:
                import traceback
                print(f"  ERROR scoring {gene}/{cond}: {e}")
                traceback.print_exc()
                continue

            pos_merged = anndata_to_per_variant_maxabs(pos_ad, pos_df)
            pos_merged["gene"] = gene
            pos_merged["variant_label"] = "pathogenic"
            pos_merged["condition"] = cond

            neg_merged = anndata_to_per_variant_maxabs(neg_ad, neg_df)
            neg_merged["gene"] = gene
            neg_merged["variant_label"] = "benign"
            neg_merged["condition"] = cond

            merged = pd.concat([pos_merged, neg_merged], ignore_index=True)
            all_per_variant.append(merged)

            # Save per-condition per-variant scores
            merged.to_csv(
                OUTPUT_DIR / f"per_variant_{gene}_{cond}.csv", index=False)

            # Per-variant aggregation (mean across 3 alts) for stats
            # Group by (pos, ref) → max-abs is averaged across alts
            pos_var = (pos_merged
                       .groupby(["pos", "ref"], as_index=False)
                       .agg({"max_abs": "mean", "mean_abs": "mean"}))
            neg_var = (neg_merged
                       .groupby(["pos", "ref"], as_index=False)
                       .agg({"max_abs": "mean", "mean_abs": "mean"}))

            n_path = len(pos_var)
            n_ben = len(neg_var)
            path_max_mean = float(pos_var["max_abs"].mean())
            ben_max_mean = float(neg_var["max_abs"].mean())
            u_max, p_max, r_max, _, _ = mannwhitney_r(
                pos_var["max_abs"], neg_var["max_abs"], alternative="greater")
            u_mean, p_mean, r_mean, _, _ = mannwhitney_r(
                pos_var["mean_abs"], neg_var["mean_abs"], alternative="greater")

            row = {
                "gene": gene,
                "condition": cond,
                "n_path_per_variant": n_path,
                "n_ben_per_variant": n_ben,
                "n_tracks": int(pos_merged["n_tracks"].dropna().iloc[0]) if pos_merged["n_tracks"].notna().any() else 0,
                "path_max_mean": path_max_mean,
                "ben_max_mean": ben_max_mean,
                "path_mean_mean": float(pos_var["mean_abs"].mean()),
                "ben_mean_mean": float(neg_var["mean_abs"].mean()),
                "u_max_abs": u_max,
                "p_max_abs": p_max,
                "r_max_abs": r_max,
                "u_mean_abs": u_mean,
                "p_mean_abs": p_mean,
                "r_mean_abs": r_mean,
                "runtime_s": pos_elapsed + neg_elapsed,
            }
            summary_rows.append(row)

            print(f"    n_path(per-var)={n_path}, n_ben(per-var)={n_ben}, "
                  f"n_tracks={row['n_tracks']}")
            print(f"    max-abs:  path={path_max_mean:.4f}  ben={ben_max_mean:.4f}  "
                  f"U={u_max:.1f}  p={p_max:.4g}  r={r_max:+.3f}  "
                  f"{'*** SIG' if p_max < 0.05 else 'ns'}")
            print(f"    mean-abs: path={row['path_mean_mean']:.4f}  "
                  f"ben={row['ben_mean_mean']:.4f}  "
                  f"U={u_mean:.1f}  p={p_mean:.4g}  r={r_mean:+.3f}  "
                  f"{'*** SIG' if p_mean < 0.05 else 'ns'}")

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTPUT_DIR / "summary_by_gene_condition.csv", index=False)
    print(f"\nSaved -> {OUTPUT_DIR / 'summary_by_gene_condition.csv'}")

    # Comparison: averaged vs tissue-filtered
    print("\n" + "=" * 70)
    print("COMPARISON: averaged vs tissue-filtered")
    print("=" * 70)
    pivot_p = summary_df.pivot_table(
        index="gene", columns="condition", values="p_max_abs")
    pivot_r = summary_df.pivot_table(
        index="gene", columns="condition", values="r_max_abs")

    print("\nMax-abs p-values:")
    print(pivot_p.to_string())
    print("\nMax-abs r (rank-biserial):")
    print(pivot_r.to_string())

    if "tissue" in pivot_p.columns and "avg" in pivot_p.columns:
        delta_p = pivot_p["tissue"] - pivot_p["avg"]
        delta_r = pivot_r["tissue"] - pivot_r["avg"]
        print("\nDelta p (tissue - avg), negative = tissue p-value SMALLER (better):")
        print(delta_p.to_string())
        print("\nDelta r (tissue - avg), positive = tissue effect LARGER:")
        print(delta_r.to_string())

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        df_plot = pd.concat(all_per_variant, ignore_index=True)
        # Per-variant aggregation for plotting
        df_plot_var = (df_plot
                       .groupby(["gene", "variant_label", "condition",
                                 "pos", "ref"], as_index=False)
                       .agg({"max_abs": "mean"}))

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, gene in zip(axes, GENES):
            sub = df_plot_var[df_plot_var["gene"] == gene]
            x_pos = []
            x_labels = []
            x_colors = []
            for i, (cond, label) in enumerate([
                ("avg", "Averaged\n(305 tracks)"),
                ("tissue", "Tissue-filtered\n(brain/muscle)"),
            ]):
                cond_data = sub[sub["condition"] == cond]
                path = cond_data[cond_data["variant_label"] == "pathogenic"]["max_abs"].dropna()
                ben = cond_data[cond_data["variant_label"] == "benign"]["max_abs"].dropna()
                # jitter
                rng = np.random.default_rng(42)
                xp = rng.normal(i - 0.18, 0.04, size=len(path))
                xb = rng.normal(i + 0.18, 0.04, size=len(ben))
                ax.scatter(xp, path, color="red", alpha=0.6, label="pathogenic" if i == 0 else None)
                ax.scatter(xb, ben, color="blue", alpha=0.6, label="benign" if i == 0 else None)
                ax.hlines(path.mean(), xmin=i - 0.32, xmax=i - 0.04,
                           colors="red", linestyles="dashed", linewidth=2)
                ax.hlines(ben.mean(), xmin=i + 0.04, xmax=i + 0.32,
                           colors="blue", linestyles="dashed", linewidth=2)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Averaged", "Tissue-filtered"])
            ax.set_xlim(-0.6, 1.6)
            ax.set_ylabel("max(|DNASE|) across tracks")
            ax.set_title(f"{gene.upper()}")
            ax.legend(loc="upper right")
            ax.grid(True, axis="y", alpha=0.3)
        fig.suptitle("Experiment 006 — DNASE: averaged vs tissue-filtered\n"
                     "(per-variant mean across 3 alt alleles)",
                     fontsize=12)
        fig.tight_layout()
        out_png = OUTPUT_DIR / "tissue_vs_averaged.png"
        fig.savefig(out_png, dpi=150)
        print(f"\nSaved -> {out_png}")
    except Exception as e:
        print(f"Plot failed: {e}")

    print(f"\nAll outputs in {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())