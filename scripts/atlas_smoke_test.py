#!/usr/bin/env python3
"""
Atlas smoke test + exploration.

Atlas is the pre-computed AlphaGenome Variant Impact (AVI) scores for all
9 billion possible SNVs. This script verifies Atlas access and shows you
what data you get back.

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/atlas_smoke_test.py
"""

from __future__ import annotations

import os
import sys

import pandas as pd

from alphagenome.atlas import atlas as atlas_mod
from alphagenome.data import genome


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    print("Creating AtlasClient...")
    atlas_client = atlas_mod.create(api_key)

    print("Listing available scorers (Atlas only supports a subset)...")
    metadata = atlas_client.scorer_metadata()
    print(f"Atlas has {len(metadata)} scorers")
    print()
    print("First 10:")
    for name, m in list(metadata.items())[:10]:
        n_tracks = len(m.track_metadata)
        print(f"  {name}: {n_tracks} tracks, signed={m.is_signed}")

    print("\nQuerying a 64-bp region on chr22 for one variant we know exists...")
    # This is the variant from quick_start.ipynb — should be in the 1000 Genomes cohort.
    test_interval = genome.Interval("chr22", 36_201_666, 36_201_730)

    print(f"Interval: {test_interval}")
    print("Calling atlas.get_interval() — chunks internally, ~30s first call...")

    scores = atlas_client.get_interval(
        interval=test_interval,
        requested_scorers=list(metadata.keys())[:3],  # First 3 scorers only
    )
    print(f"Got AnnData objects for {len(scores)} scorers:")
    for scorer_name, adata in scores.items():
        print(f"\n--- {scorer_name} ---")
        print(f"Shape: {adata.X.shape}  → (n_variants, n_tracks)")
        print(f"Layers: {list(adata.layers.keys()) if adata.layers else 'none'}")
        if adata.obs is not None and len(adata.obs) > 0:
            print(f"\n.obs (variant metadata) head:")
            print(adata.obs.head().to_string())
        if adata.var is not None and len(adata.var) > 0:
            print(f"\n.var (track metadata) head:")
            print(adata.var.head().to_string())

    print("\nDONE")
    print("If you reached here, Atlas access works and you understand the data shape.")
    print("Next: explore what scorers are useful for your disease of interest.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
