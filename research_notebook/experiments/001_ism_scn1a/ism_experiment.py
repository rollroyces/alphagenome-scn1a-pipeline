#!/usr/bin/env python3
"""
Experiment 001: ISM (in-silico mutagenesis) on SCN1A variants.

Goal: discover what sequence features AlphaGenome uses to score pathogenic
splice-disrupting variants differently from benign controls.

Hypothesis:
- Pathogenic variants' ISM profiles show concentrated sensitivity at canonical
  splice sites (positions ±1, ±2 of intron-exon boundaries) AND at cryptic
  splice activators deeper in introns.
- Benign variants' ISM profiles show diffuse, weaker sensitivity.
- This difference should be visible as a concentration metric (e.g., the
  fraction of total ISM mass within ±5 bp of the variant position).

If we see this pattern: confirms AlphaGenome has learned canonical biology.
If we see something different: potential discovery.

Approach:
- Pick 10 pathogenic splicing variants from our benchmark (highest scores)
- Pick 10 benign intronic controls (lowest scores, matched distance from exons)
- Run ISM with a 64 bp window centered on each variant
- Compute ISM concentration metrics, plot heatmaps

Runtime: ~15-30 minutes for 20 variants × 64 bp × 3 alts = 3,840 variant scores
Output: research_notebook/experiments/001_ism_scn1a/
"""

from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# Configuration
ISM_WINDOW = 64  # bp on each side of the variant
CONTEXT = dna_client.SEQUENCE_LENGTH_16KB  # 16,384 bp context
N_PATHOGENIC = 10
N_BENIGN = 10
OUTPUT_DIR = Path("research_notebook/experiments/001_ism_scn1a")

# Smoke test mode — only run 1+1 variants for fast iteration
import os
SMOKE_TEST = os.environ.get("ISM_SMOKE") == "1"
if SMOKE_TEST:
    N_PATHOGENIC = 1
    N_BENIGN = 1


def load_benchmark_variants():
    """Pick top-N pathogenic and bottom-N benign intronic variants from SCN1A benchmark (SNVs only)."""
    bench = pd.read_csv("outputs/benchmark_scn1a_live_api_raw.csv")
    # Filter to SNVs only (ref and alt must each be 1 bp) — AlphaGenome variant scoring is SNV-focused
    bench = bench[(bench["ref"].str.len() == 1) & (bench["alt"].str.len() == 1)].copy()
    pos = bench[bench["clnsig_category"] == "pathogenic"].sort_values("SPLICE_SITES_score", ascending=False)
    neg = bench[bench["clnsig_category"] == "benign"].sort_values("SPLICE_SITES_score", ascending=True)
    return pos.head(N_PATHOGENIC).reset_index(drop=True), neg.head(N_BENIGN).reset_index(drop=True)


def run_ism(dna_model, chrom: str, pos: int, ref: str, alt: str, variant_label: str):
    """Run ISM on a single variant position.

    Returns: ISM matrix (2*ISM_WINDOW, 4) with score delta per (position, base).
    """
    # Center the sequence on the variant
    interval = genome.Interval(chrom, pos - CONTEXT // 2, pos + CONTEXT // 2)
    # Window to ISM: ISM_WINDOW bp centered on variant
    ism_interval = genome.Interval(chrom, pos - ISM_WINDOW, pos + ISM_WINDOW)

    # Scorer: SPLICE_SITES (matches what we benchmarked)
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"]

    variant_scores = dna_model.score_ism_variants(
        interval=interval,
        ism_interval=ism_interval,
        variant_scorers=[scorer],
    )

    # variant_scores is list[list[AnnData]] — outer list = variants, inner list = scorers
    # We use one scorer, so flatten to scores per variant
    ism_values = []
    for variant_scores_list in variant_scores:
        # variant_scores_list is a list of AnnData (one per scorer)
        # Sum across all scorers' tracks
        total = sum(float(ann.X.sum()) for ann in variant_scores_list)
        ism_values.append(total)

    return ism_values, ism_interval


def build_ism_matrix(ism_values, variants):
    """Construct ISM matrix from variant scores."""
    # variants: list of genome.Variant objects
    # ism_values: list of floats, same length
    # Returns: matrix of shape (n_positions, 4) with score for each (pos, alt_base)
    bases = ['A', 'C', 'G', 'T']
    base_idx = {b: i for i, b in enumerate(bases)}
    n_positions = ISM_WINDOW * 2

    matrix = np.zeros((n_positions, 4), dtype=np.float32)
    for v, s in zip(variants, ism_values):
        # variant.position is 1-based
        rel_pos = v.position - (min(va.position for va in variants))
        if 0 <= rel_pos < n_positions:
            alt = v.alternate_bases
            if alt in base_idx:
                matrix[rel_pos, base_idx[alt]] = s
    return matrix


def compute_concentration_metrics(ism_matrix):
    """Compute how concentrated the ISM signal is near the variant (center)."""
    n_positions = ism_matrix.shape[0]
    center = n_positions // 2

    # Compute magnitude (sum of |effects| across bases at each position)
    magnitude = np.abs(ism_matrix).sum(axis=1)

    # Total effect
    total = magnitude.sum()
    if total == 0:
        return {"total": 0, "fraction_within_5bp": 0, "fraction_within_15bp": 0,
                "concentration_5bp": 0, "max_position": -1, "max_value": 0}

    # Concentration near center
    within_5bp = magnitude[center - 5: center + 6].sum()  # ±5 bp = 11 positions
    within_15bp = magnitude[center - 15: center + 16].sum()  # ±15 bp = 31 positions

    # Position of max effect
    max_pos = int(np.argmax(magnitude))
    max_val = float(magnitude[max_pos])

    return {
        "total": float(total),
        "fraction_within_5bp": float(within_5bp / total),
        "fraction_within_15bp": float(within_15bp / total),
        "concentration_5bp": float(within_5bp / (11 * magnitude.mean())) if magnitude.mean() > 0 else 0,
        "max_position_relative_to_variant": int(max_pos - center),
        "max_value": max_val,
    }


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

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    all_results = []

    for label, df_subset, n_total in [
        ("pathogenic", pos_variants, N_PATHOGENIC),
        ("benign", neg_variants, N_BENIGN),
    ]:
        print(f"\n--- Running ISM on {n_total} {label} variants ---")

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
                ism_values, ism_interval = run_ism(dna_model, chrom, pos, ref, alt, label)
                # Convert to numpy matrix
                ism_array = np.array(ism_values)
                # Each variant has 3 alternative base scores → 3 values per position
                n_positions = ISM_WINDOW * 2
                # Reshape into (n_positions, 3) matrix (alt bases only)
                if len(ism_array) != n_positions * 3:
                    raise ValueError(f"Got {len(ism_array)} scores, expected {n_positions * 3}")
                matrix = ism_array.reshape(n_positions, 3)

                # Construct full (n_positions, 4) matrix with ref base = 0
                bases = ['A', 'C', 'G', 'T']
                alts = [b for b in bases if b != ref]
                full_matrix = np.zeros((n_positions, 4), dtype=np.float32)
                for j, b in enumerate(alts):
                    full_matrix[:, bases.index(b)] = matrix[:, j]

                metrics = compute_concentration_metrics(full_matrix)
                metrics["variant_label"] = label
                metrics["chrom"] = chrom
                metrics["pos"] = pos
                metrics["ref"] = ref
                metrics["alt"] = alt
                metrics["rank"] = i
                metrics["runtime_s"] = time.time() - t0
                metrics["success"] = True
                all_results.append(metrics)

                if i % 5 == 0:
                    print(f"  [{i+1}/{n_total}] {metrics['fraction_within_5bp']:.3f} in ±5bp, "
                          f"max@{metrics['max_position_relative_to_variant']:+d}bp "
                          f"({metrics['runtime_s']:.0f}s)", flush=True)

                # Save ISM matrix for later visualization
                np.save(OUTPUT_DIR / f"ism_{label}_{i:02d}_pos{pos}.npy", full_matrix)

            except Exception as e:
                print(f"  [{i+1}/{n_total}] FAILED: {type(e).__name__}: {str(e)[:80]}", flush=True)
                all_results.append({
                    "variant_label": label, "rank": i, "success": False,
                    "error": f"{type(e).__name__}: {str(e)[:100]}",
                })

    # Save results
    df = pd.DataFrame(all_results)
    df.to_csv(OUTPUT_DIR / "metrics.csv", index=False)
    print(f"\nSaved → {OUTPUT_DIR / 'metrics.csv'}")

    # Summary statistics
    print("\n" + "=" * 60)
    print("ISM CONCENTRATION RESULTS")
    print("=" * 60)
    for label in ["pathogenic", "benign"]:
        sub = df[(df["variant_label"] == label) & (df["success"] == True)]
        if len(sub) == 0:
            continue
        print(f"\n{label.upper()} (n={len(sub)}):")
        print(f"  Fraction of effect within ±5bp of variant: "
              f"{sub['fraction_within_5bp'].mean():.3f} ± {sub['fraction_within_5bp'].std():.3f}")
        print(f"  Fraction within ±15bp: "
              f"{sub['fraction_within_15bp'].mean():.3f} ± {sub['fraction_within_15bp'].std():.3f}")
        print(f"  Total ISM magnitude: "
              f"{sub['total'].mean():.3f} ± {sub['total'].std():.3f}")

    # Statistical test: are pathogenic variants more concentrated?
    pos_sub = df[(df["variant_label"] == "pathogenic") & (df["success"] == True)]
    neg_sub = df[(df["variant_label"] == "benign") & (df["success"] == True)]
    if len(pos_sub) >= 3 and len(neg_sub) >= 3:
        from scipy.stats import mannwhitneyu
        stat, p = mannwhitneyu(
            pos_sub["fraction_within_5bp"].dropna(),
            neg_sub["fraction_within_5bp"].dropna(),
            alternative="two-sided"
        )
        print(f"\nMann-Whitney U test (pathogenic vs benign ±5bp concentration):")
        print(f"  U = {stat}, p = {p:.4f}")
        if p < 0.05:
            print(f"  → SIGNIFICANT: pathogenic variants show different ISM concentration")
        else:
            print(f"  → Not significant at α=0.05")

    print(f"\nAll outputs in {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
