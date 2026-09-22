#!/usr/bin/env python3
"""
In-silico mutagenesis (ISM) example from quick_start.ipynb.

What ISM does: given a region, systematically mutate every base to every other
base and observe the prediction change. This identifies which positions are
functionally important.

This is the operation that powers most "non-coding variant prioritization"
use cases. Understanding ISM is critical for the rare-disease work.

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/ism_example.py
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    dna_model = dna_client.create(api_key)

    print("Setting up ISM region: chr20:3,753,000–3,753,400 (400 bp), 16 Kb context")
    sequence_interval = genome.Interval("chr20", 3_753_000, 3_753_400)
    sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)

    # Mutate every base in the central 256 bp of the sequence_interval.
    ism_interval = sequence_interval.resize(256)

    print(f"ISM interval (where we mutate every base): {ism_interval}")
    print(f"Number of variants generated: 256 bp × 3 alt alleles = 768 mutations")

    # Score each mutation's effect on a 501 bp window of DNase centered on the mutation.
    dnase_variant_scorer = variant_scorers.CenterMaskScorer(
        requested_output=dna_client.OutputType.DNASE,
        width=501,
        aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
    )
    print(f"Scorer: {dnase_variant_scorer}")

    print("\nCalling API for 768 variants — this may take 30–60 seconds...")
    variant_scores = dna_model.score_ism_variants(
        interval=sequence_interval,
        ism_interval=ism_interval,
        variant_scorers=[dnase_variant_scorer],
    )
    print(f"Got {len(variant_scores)} AnnData objects")

    adata = variant_scores[0]
    print(f"AnnData shape: {adata.X.shape}")
    print(f"  → 768 variants × N tracks (tissue/cell types)")

    # Extract K562 score (a well-studied cell line) for visualization.
    def extract_k562(a):
        values = a.X[:, a.var["ontology_curie"] == "EFO:0002067"]
        assert values.size == 1
        return values.flatten()[0]

    # Build a (256 bp × 3 alts) matrix of scores
    ism_matrix = np.full((ism_interval.width, 3), np.nan)
    ref_seq = str(dna_model.reference_genome.fetch(ism_interval)).upper()
    print(f"Reference sequence ({len(ref_seq)} bp): {ref_seq}")

    for ad in variant_scores:
        row = extract_k562(ad)
        var_pos = ad.obs["position"].iloc[0] - ism_interval.start
        alt = ad.obs["alt"].iloc[0]
        alt_idx = "ACGT".index(alt)
        if 0 <= var_pos < ism_interval.width:
            ism_matrix[var_pos, alt_idx] = row

    # Plot: heatmap of K562 DNase effect for each (position, alt) combination.
    fig, ax = plt.subplots(figsize=(14, 4))
    vmax = np.nanpercentile(np.abs(ism_matrix), 95)
    im = ax.imshow(
        ism_matrix,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_xlabel("Alt allele (A, C, G, T)")
    ax.set_ylabel(f"Position in {ism_interval}")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["A", "C", "G"])
    ax.set_title("K562 DNase ISM effect — chr20:3,753,072–3,753,328")
    plt.colorbar(im, ax=ax, label="Δ DNase prediction")
    plt.tight_layout()

    out_path = "figures/ism_k562_heatmap.png"
    os.makedirs("figures", exist_ok=True)
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved ISM heatmap → {out_path}")

    # Highlight the top 5 most impactful positions
    impact = np.nanmax(np.abs(ism_matrix), axis=1)
    top5 = np.argsort(impact)[-5:][::-1]
    print("\nTop 5 most impactful positions in the ISM region:")
    for pos in top5:
        pos_actual = ism_interval.start + pos
        ref = ref_seq[pos]
        alts = "ACGT"
        impacts = []
        for alt_idx in range(3):
            alts_ = alts[alt_idx] if alt_idx < 2 else alts[alt_idx + (1 if alt_idx >= 2 else 0)]
            # simpler:
        for alt_idx in range(3):
            alt_letter = alts[alt_idx]
            if alt_letter == ref:
                continue
            score = ism_matrix[pos, alt_idx]
            if not np.isnan(score):
                impacts.append(f"{ref}→{alt_letter}: {score:+.4f}")
        print(f"  pos {pos_actual:,} ({ref}): {' | '.join(impacts)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
