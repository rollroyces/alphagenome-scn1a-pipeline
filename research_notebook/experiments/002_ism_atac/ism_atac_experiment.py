#!/usr/bin/env python3
"""
Experiment 002: ISM on chromatin-level modalities (ATAC + DNase) for SCN1A.

Goal: discover whether ATAC and DNase show a DIFFERENT concentration pattern
than splicing (Experiment 001) for pathogenic vs benign SCN1A variants.

Hypothesis:
- Experiment 001 (splicing) found pathogenic and benign variants have SIMILAR
  spatial concentration; only magnitude differs.
- Chromatin-level features (ATAC, DNase) might show a different pattern,
  because AlphaGenome has to integrate broader regulatory context for them.

Approach:
- Same 10 pathogenic + 10 benign SCN1A intronic SNVs as Experiment 001.
- Run ISM with the recommended ATAC CenterMaskScorer (width 501, DIFF_LOG2_SUM).
- ALSO run with DNase CenterMaskScorer for direct comparison.
- ISM window: 64 bp on each side of the variant (128 bp total), 16,384 bp context.
- Compute same metrics: fraction within +/-5 bp, +/-15 bp, total magnitude.

Runtime: ~5-10 minutes per modality for 20 variants x 128 bp x 3 alts
= 7,680 variant scores x 2 scorers = ~30-50 minutes total.

Output: research_notebook/experiments/002_ism_atac/

Usage:
    export ISM_SMOKE=1 && bash scripts/_run_with_key.sh research_notebook/experiments/002_ism_atac/ism_atac_experiment.py
    unset ISM_SMOKE && bash scripts/_run_with_key.sh research_notebook/experiments/002_ism_atac/ism_atac_experiment.py
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
N_PATHOGENIC = 10
N_BENIGN = 10
OUTPUT_DIR = Path("research_notebook/experiments/002_ism_atac")

# Smoke test mode - only run 1+1 variants for fast iteration
SMOKE_TEST = os.environ.get("ISM_SMOKE") == "1"
if SMOKE_TEST:
    N_PATHOGENIC = 1
    N_BENIGN = 1

# Modalities to compare
MODALITIES = [
    ("ATAC", variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"]),
    ("DNASE", variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]),
]


def load_benchmark_variants():
    """Pick top-N pathogenic and bottom-N benign intronic SNVs from SCN1A benchmark."""
    bench = pd.read_csv("outputs/benchmark_scn1a_live_api_raw.csv")
    # Filter to SNVs only (ref and alt must each be 1 bp) - AlphaGenome variant scoring is SNV-focused
    bench = bench[(bench["ref"].str.len() == 1) & (bench["alt"].str.len() == 1)].copy()
    pos = bench[bench["clnsig_category"] == "pathogenic"].sort_values("SPLICE_SITES_score", ascending=False)
    neg = bench[bench["clnsig_category"] == "benign"].sort_values("SPLICE_SITES_score", ascending=True)
    return pos.head(N_PATHOGENIC).reset_index(drop=True), neg.head(N_BENIGN).reset_index(drop=True)


def normalize_variant(v):
    """Return (chrom, pos, ref, alt) from either a Variant object or a 'chr:pos:ref>alt' string."""
    if isinstance(v, str):
        m = re.match(r"(.+):(\d+):([ACGT])>([ACGT])", v)
        if not m:
            raise ValueError(f"Cannot parse variant string: {v}")
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    # Variant object - str(v) gives "chr2:166054637:A>C" format
    return str(v.chromosome), int(v.position), str(v.reference_bases), str(v.alternate_bases)


def run_ism(dna_model, chrom: str, pos: int, ref: str, alt: str, scorer):
    """Run ISM on a single variant position.

    Returns:
        List of (variant_str_or_obj, total_score) tuples, where total_score is
        the sum across all output tracks for that variant.
    """
    interval = genome.Interval(chrom, pos - CONTEXT // 2, pos + CONTEXT // 2)
    ism_interval = genome.Interval(chrom, pos - ISM_WINDOW, pos + ISM_WINDOW)

    variant_scores = dna_model.score_ism_variants(
        interval=interval,
        ism_interval=ism_interval,
        variant_scorers=[scorer],
    )

    # variant_scores is list[list[AnnData]] - outer list = variants, inner list = scorers
    results = []
    for variant_scores_list in variant_scores:
        for ann in variant_scores_list:
            variant_obj = ann.uns.get("variant", "")
            total = float(np.asarray(ann.X).sum())
            results.append((variant_obj, total))
    return results


def build_ism_matrix(variant_results, center_pos, ref_base):
    """Construct ISM matrix (n_positions, 4) from variant results.

    For each position, we have 3 alt variants. Ref base column = 0.
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
    """Compute how concentrated the ISM signal is near the variant (center)."""
    n_positions = ism_matrix.shape[0]
    center = n_positions // 2

    # Magnitude = sum of |effects| across bases at each position
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

    within_5bp = float(magnitude[center - 5: center + 6].sum())  # +/-5 bp = 11 positions
    within_15bp = float(magnitude[center - 15: center + 16].sum())  # +/-15 bp = 31 positions

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


def run_modality(dna_model, modality_name, scorer, pos_variants, neg_variants, all_results):
    """Run ISM for one modality across all variants."""
    print(f"\n{'='*60}")
    print(f"MODALITY: {modality_name} ({scorer})")
    print(f"{'='*60}")

    for label, df_subset, n_total in [
        ("pathogenic", pos_variants, N_PATHOGENIC),
        ("benign", neg_variants, N_BENIGN),
    ]:
        print(f"\n--- Running {modality_name} ISM on {n_total} {label} variants ---")

        for i in range(n_total):
            row = df_subset.iloc[i]
            chrom = str(row["chrom"])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            pos = int(row["pos"])
            ref = str(row["ref"])
            alt = str(row["alt"])

            t0 = time.time()
            try:
                results = run_ism(dna_model, chrom, pos, ref, alt, scorer)
                # Build (n_positions, 4) matrix once per variant
                full_matrix = build_ism_matrix(results, pos, ref)

                # Validate shape
                expected = ISM_WINDOW * 2
                if full_matrix.shape != (expected, 4):
                    raise ValueError(f"Got matrix shape {full_matrix.shape}, expected ({expected}, 4)")

                metrics = compute_concentration_metrics(full_matrix)
                metrics["variant_label"] = label
                metrics["modality"] = modality_name
                metrics["chrom"] = chrom
                metrics["pos"] = pos
                metrics["ref"] = ref
                metrics["alt"] = alt
                metrics["rank"] = i
                metrics["runtime_s"] = time.time() - t0
                metrics["success"] = True
                # Single row per (variant, modality), ref/alt preserved for traceability
                all_results.append(metrics)

                print(f"  [{i+1}/{n_total}] +/-5bp frac={metrics['fraction_within_5bp']:.3f}, "
                      f"max@{metrics['max_position_relative_to_variant']:+d}bp, "
                      f"total={metrics['total']:.2f} ({metrics['runtime_s']:.0f}s)",
                      flush=True)

                # Save ISM matrix (one file per (variant, modality))
                np.save(OUTPUT_DIR / f"ism_{modality_name.lower()}_{label}_{i:02d}_pos{pos}.npy", full_matrix)

            except Exception as e:
                print(f"  [{i+1}/{n_total}] FAILED: {type(e).__name__}: {str(e)[:120]}", flush=True)
                all_results.append({
                    "variant_label": label,
                    "modality": modality_name,
                    "rank": i,
                    "success": False,
                    "error": f"{type(e).__name__}: {str(e)[:200]}",
                })


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading benchmark variants...")
    pos_variants, neg_variants = load_benchmark_variants()
    print(f"  Top {N_PATHOGENIC} pathogenic: positions {list(pos_variants['pos'].head(5))}...")
    print(f"  Bottom {N_BENIGN} benign: positions {list(neg_variants['pos'].head(5))}...")
    if SMOKE_TEST:
        print("  [SMOKE TEST MODE - 1+1 variants only]")

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    all_results = []

    # Run both modalities
    for modality_name, scorer in MODALITIES:
        run_modality(dna_model, modality_name, scorer, pos_variants, neg_variants, all_results)

    # Save results
    df = pd.DataFrame(all_results)
    df.to_csv(OUTPUT_DIR / "metrics.csv", index=False)
    print(f"\nSaved -> {OUTPUT_DIR / 'metrics.csv'}")

    # Summary by modality and label
    print("\n" + "=" * 70)
    print("ISM CONCENTRATION RESULTS BY MODALITY")
    print("=" * 70)

    from scipy.stats import mannwhitneyu

    for modality_name, _ in MODALITIES:
        print(f"\n--- {modality_name} ---")
        for label in ["pathogenic", "benign"]:
            sub = df[(df["variant_label"] == label) &
                     (df["modality"] == modality_name) &
                     (df["success"] == True)]
            if len(sub) == 0:
                print(f"  {label}: no successful variants")
                continue
            print(f"  {label.upper()} (n={len(sub)}):")
            print(f"    Fraction +/-5bp:     {sub['fraction_within_5bp'].mean():.3f} +/- {sub['fraction_within_5bp'].std():.3f}")
            print(f"    Fraction +/-15bp:    {sub['fraction_within_15bp'].mean():.3f} +/- {sub['fraction_within_15bp'].std():.3f}")
            print(f"    Total magnitude:     {sub['total'].mean():.2f} +/- {sub['total'].std():.2f}")

        # Statistical test
        pos_sub = df[(df["variant_label"] == "pathogenic") &
                     (df["modality"] == modality_name) &
                     (df["success"] == True)]
        neg_sub = df[(df["variant_label"] == "benign") &
                     (df["modality"] == modality_name) &
                     (df["success"] == True)]
        if len(pos_sub) >= 3 and len(neg_sub) >= 3:
            stat, p = mannwhitneyu(
                pos_sub["fraction_within_5bp"].dropna(),
                neg_sub["fraction_within_5bp"].dropna(),
                alternative="two-sided",
            )
            print(f"  Mann-Whitney U (+/-5bp): U={stat}, p={p:.4f}", end="")
            if p < 0.05:
                print(f" -> SIGNIFICANT")
            else:
                print(f" -> not significant")

            stat_mag, p_mag = mannwhitneyu(
                pos_sub["total"].dropna(),
                neg_sub["total"].dropna(),
                alternative="two-sided",
            )
            print(f"  Mann-Whitney U (magnitude): U={stat_mag}, p={p_mag:.4f}", end="")
            if p_mag < 0.05:
                print(f" -> SIGNIFICANT")
            else:
                print(f" -> not significant")

    print(f"\nAll outputs in {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
