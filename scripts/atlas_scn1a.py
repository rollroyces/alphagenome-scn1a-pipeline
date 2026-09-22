#!/usr/bin/env python3
"""
Query AlphaGenome Atlas for ALL SCN1A variants — splicing scores only.

Atlas is the pre-computed AlphaGenome Variant Impact (AVI) scores for every
possible SNV in the human genome. This is the cheap path to genome-wide
scoring (no live API quota consumed).

We query just the SCN1A locus, filtered to splicing scorers. Output is a
DataFrame with one row per (variant, allele) pair, with predicted splicing
impact scores.

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/atlas_scn1a.py

Output:
    data/atlas_scn1a_splicing.h5ad — AnnData with all splicing scores
    outputs/atlas_scn1a_summary.csv — flat CSV for downstream analysis
"""

from __future__ import annotations

import os
import sys

import anndata
import numpy as np
import pandas as pd

from alphagenome.atlas import atlas
from alphagenome.data import genome


# SCN1A gene interval (hg38, minus strand)
SCN1A_CHROM = "chr2"
SCN1A_START = 165_984_640
SCN1A_END = 166_182_806
SCN1A_FLANK = 50_000  # regulatory variants can extend into flanking regions


# Splicing scorers — these are the ones relevant for SCN1A noncoding variants
# NOTE: SPLICE_JUNCTIONS_ACTIVE is not in Atlas (verified empirically)
SPLICING_SCORERS = [
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
]


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    # Monkey-patch Atlas gRPC channel to allow large response messages.
    # The default gRPC max is 4 MB; some Atlas responses for ~32 bp chunks ×
    # 367 tracks × multiple scorers exceed this.
    import grpc
    _orig_secure_channel = grpc.secure_channel

    def _patched_secure_channel(address, credentials, options=()):
        new_options = list(options) + [
            ("grpc.max_receive_message_length", -1),
            ("grpc.max_send_message_length", -1),
        ]
        return _orig_secure_channel(address, credentials, options=tuple(new_options))

    grpc.secure_channel = _patched_secure_channel
    print("Patched grpc.secure_channel for unlimited message size")

    start = max(0, SCN1A_START - SCN1A_FLANK)
    end = SCN1A_END + SCN1A_FLANK
    interval = genome.Interval(SCN1A_CHROM, start, end)
    print(f"Querying Atlas for region {interval} (SCN1A ± {SCN1A_FLANK:,} bp)")
    print(f"Total span: {(end - start):,} bp → ~{(end - start) * 3:,} variants")

    print("Creating AtlasClient...")
    atlas_client = atlas.create(api_key)

    # First check which scorers Atlas actually has
    print("\nDiscovering Atlas scorers...")
    metadata = atlas_client.scorer_metadata()
    available = set(metadata.keys())
    print(f"Atlas has {len(available)} scorers total")

    # Filter to the splicing scorers we want
    splicing_to_query = [s for s in SPLICING_SCORERS if s in available]
    missing = [s for s in SPLICING_SCORERS if s not in available]
    if missing:
        print(f"⚠ Scorers not available in Atlas: {missing}")
    print(f"Will query {len(splicing_to_query)} splicing scorers: {splicing_to_query}")

    # Show what each requested scorer contains
    for scorer_name in splicing_to_query:
        m = metadata[scorer_name]
        n_tracks = len(m.track_metadata)
        print(f"  {scorer_name}: {n_tracks} tracks, signed={m.is_signed}")
        if n_tracks > 0:
            print(f"    First 3 tracks: {m.track_metadata.index[:3].tolist()}")

    # Query
    print(f"\nCalling atlas.query_interval()...")
    print(f"This will take ~1–5 minutes (depends on chunking + network)")
    print(f"Estimated chunks: {(end - start) // 32:,}")

    scores_by_scorer = atlas_client.query_interval(
        interval=interval,
        requested_scorers=splicing_to_query,
        # Atlas default is 10. Higher values cause thread pool / gRPC channel
        # exhaustion — semaphore leaks at ~50% progress. 8 is safer.
        max_workers=8,
    )

    print(f"\nGot scores for {len(scores_by_scorer)} scorers:")
    for scorer_name, adata in scores_by_scorer.items():
        print(f"  {scorer_name}: shape {adata.X.shape}, "
              f"obs={adata.obs.shape if adata.obs is not None else 'None'}, "
              f"var={adata.var.shape if adata.var is not None else 'None'}")

    # Save full AnnData objects
    os.makedirs("data", exist_ok=True)
    for scorer_name, adata in scores_by_scorer.items():
        out_path = f"data/atlas_scn1a_{scorer_name.lower()}.h5ad"
        adata.write_h5ad(out_path)
        print(f"Saved → {out_path}")

    # Build a flat CSV summary
    # Each row = one (variant_position, alt_allele) combination
    # Columns: pos, ref, alt, scorer1_score, scorer2_score, ...
    print("\nBuilding flat summary CSV...")

    # Use SPLICE_JUNCTIONS as the primary key (it has per-variant rows)
    primary_scorer = "SPLICE_JUNCTIONS" if "SPLICE_JUNCTIONS" in scores_by_scorer else list(scores_by_scorer.keys())[0]
    primary = scores_by_scorer[primary_scorer]

    rows = []
    if primary.obs is not None and "variant" in primary.obs.columns:
        for idx in range(primary.shape[0]):
            variant = primary.obs.iloc[idx]["variant"]
            row = {
                "chrom": variant.chromosome,
                "pos": variant.position,
                "ref": variant.reference_bases,
                "alt": variant.alternate_bases,
            }
            # Add score from each scorer
            for scorer_name, adata in scores_by_scorer.items():
                # Find the matching row (same variant) in this scorer
                if adata.obs is not None and "variant" in adata.obs.columns:
                    # Match by (pos, ref, alt)
                    matches = adata.obs[
                        (adata.obs["variant"].apply(lambda v: v.position) == variant.position)
                        & (adata.obs["variant"].apply(lambda v: v.reference_bases) == variant.reference_bases)
                        & (adata.obs["variant"].apply(lambda v: v.alternate_bases) == variant.alternate_bases)
                    ]
                    if len(matches) > 0:
                        # Score is the X value at the matching row, summed across tracks
                        match_idx = matches.index[0]
                        scorer_score = adata.X[adata.obs.index.get_loc(match_idx)].sum()
                        row[f"{scorer_name}_score"] = float(scorer_score)
                        # Calibrated quantile if available
                        if adata.layers is not None and "quantiles" in adata.layers:
                            q_score = adata.layers["quantiles"][adata.obs.index.get_loc(match_idx)].sum()
                            row[f"{scorer_name}_quantile"] = float(q_score)
            rows.append(row)

    df = pd.DataFrame(rows)
    print(f"Built summary: {df.shape}")

    out_csv = "outputs/atlas_scn1a_summary.csv"
    os.makedirs("outputs", exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Saved → {out_csv}")

    # Quick summary stats
    print("\n=== Score distribution per scorer ===")
    for col in df.columns:
        if col.endswith("_score"):
            vals = df[col].dropna()
            if len(vals) > 0:
                print(f"  {col}: mean={vals.mean():.4f}, std={vals.std():.4f}, "
                      f"min={vals.min():.4f}, max={vals.max():.4f}, "
                      f"top 1% threshold={vals.quantile(0.99):.4f}")

    print("\nDONE")
    print(f"Next: cross-reference with ClinVar in scripts/benchmark_scn1a.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
