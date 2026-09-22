#!/usr/bin/env python3
"""
Query AlphaGenome Atlas for SCN1A — split into chunks to avoid SDK exhaustion.

The full locus query (298 Kb) consistently fails at ~47% with a gRPC semaphore
leak in the Atlas SDK. Workaround: split into ~30 Kb sub-intervals, query each
sequentially, save after each so partial progress survives crashes.

Output:
- data/atlas_scn1a_<scorer>.h5ad per sub-interval
- data/atlas_scn1a_combined_<scorer>.h5ad merged at end
- outputs/atlas_scn1a_summary.csv flat CSV for downstream analysis

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/atlas_scn1a_chunked.py
"""

from __future__ import annotations

import os
import sys
import time

import anndata as ad
import pandas as pd

from alphagenome.atlas import atlas
from alphagenome.data import genome


SCN1A_CHROM = "chr2"
SCN1A_START = 165_984_640
SCN1A_END = 166_182_806
SCN1A_FLANK = 50_000
CHUNK_SIZE = 30_000  # 30 Kb sub-intervals — empirically safe

SPLICING_SCORERS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    # Monkey-patch grpc.secure_channel for large responses
    import grpc

    _orig = grpc.secure_channel

    def _patched(address, credentials, options=()):
        return _orig(
            address,
            credentials,
            options=tuple(list(options) + [
                ("grpc.max_receive_message_length", -1),
                ("grpc.max_send_message_length", -1),
            ]),
        )

    grpc.secure_channel = _patched

    start = max(0, SCN1A_START - SCN1A_FLANK)
    end = SCN1A_END + SCN1A_FLANK
    total_span = end - start
    print(f"Total span: {total_span:,} bp")
    print(f"Chunk size: {CHUNK_SIZE:,} bp → {((total_span - 1) // CHUNK_SIZE) + 1} chunks")

    atlas_client = atlas.create(api_key)
    metadata = atlas_client.scorer_metadata()
    available = set(metadata.keys())
    splicing_to_query = [s for s in SPLICING_SCORERS if s in available]
    print(f"Splicing scorers available: {splicing_to_query}")

    os.makedirs("data/atlas_chunks", exist_ok=True)

    # Track all variant records we've seen
    all_records: dict[str, list[ad.AnnData]] = {s: [] for s in splicing_to_query}
    processed = 0
    t0 = time.time()

    chunk_start = start
    chunk_idx = 0
    while chunk_start < end:
        chunk_end = min(chunk_start + CHUNK_SIZE, end)
        interval = genome.Interval(SCN1A_CHROM, chunk_start, chunk_end)
        chunk_idx += 1
        n_chunks = ((total_span - 1) // CHUNK_SIZE) + 1
        print(f"\n[{chunk_idx}/{n_chunks}] Querying {interval} ({interval.width:,} bp)...")

        try:
            scores_by_scorer = atlas_client.query_interval(
                interval=interval,
                requested_scorers=splicing_to_query,
                max_workers=4,
            )
            for scorer_name, ann_data in scores_by_scorer.items():
                all_records[scorer_name].append(ann_data)
                # Save chunk for resumability
                chunk_path = f"data/atlas_chunks/{scorer_name}_chunk_{chunk_idx:03d}.h5ad"
                ann_data.write_h5ad(chunk_path)
            processed += interval.width
            elapsed = time.time() - t0
            rate = processed / elapsed if elapsed > 0 else 0
            remaining = (total_span - processed) / rate if rate > 0 else 0
            print(f"  ✓ Saved {len(scores_by_scorer)} scorers")
            print(f"  Progress: {processed:,}/{total_span:,} bp ({100*processed/total_span:.1f}%) — "
                  f"{rate:.0f} bp/s — ETA {remaining:.0f}s")
        except Exception as e:
            print(f"  ✗ Failed: {type(e).__name__}: {e}")
            print(f"  Skipping chunk, continuing...")
            # Clean shutdown to avoid semaphore leak carrying forward
            try:
                del atlas_client
            except Exception:
                pass
            atlas_client = atlas.create(api_key)
            metadata = atlas_client.scorer_metadata()

        chunk_start = chunk_end

    print("\n" + "=" * 70)
    print("All chunks processed. Merging...")
    print("=" * 70)

    combined: dict[str, ad.AnnData] = {}
    for scorer_name, adatas in all_records.items():
        if not adatas:
            print(f"  {scorer_name}: no data, skipping")
            continue
        merged = ad.concat(adatas, axis=0, label="chunk_idx", index_unique=None)
        combined[scorer_name] = merged
        print(f"  {scorer_name}: merged {len(adatas)} chunks → {merged.shape}")

        out_path = f"data/atlas_scn1a_combined_{scorer_name.lower()}.h5ad"
        merged.write_h5ad(out_path)
        print(f"    Saved → {out_path}")

    if not combined:
        print("ERROR: no data collected")
        return 2

    # Build flat CSV summary from the combined SPLICE_JUNCTIONS data
    print("\nBuilding flat CSV summary...")
    primary_scorer = "SPLICE_JUNCTIONS" if "SPLICE_JUNCTIONS" in combined else list(combined.keys())[0]
    primary = combined[primary_scorer]

    # Build a per-variant dict keyed by (pos, ref, alt)
    if primary.obs is not None and "variant" in primary.obs.columns:
        rows = []
        for idx in range(primary.shape[0]):
            variant = primary.obs.iloc[idx]["variant"]
            row = {
                "chrom": variant.chromosome,
                "pos": variant.position,
                "ref": variant.reference_bases,
                "alt": variant.alternate_bases,
            }
            for scorer_name, adata in combined.items():
                if adata.obs is not None and "variant" in adata.obs.columns:
                    matches = adata.obs[
                        (adata.obs["variant"].apply(lambda v: v.position) == variant.position)
                        & (adata.obs["variant"].apply(lambda v: v.reference_bases) == variant.reference_bases)
                        & (adata.obs["variant"].apply(lambda v: v.alternate_bases) == variant.alternate_bases)
                    ]
                    if len(matches) > 0:
                        match_idx = matches.index[0]
                        adata_idx = adata.obs.index.get_loc(match_idx)
                        row[f"{scorer_name}_score"] = float(adata.X[adata_idx].sum())
                        if adata.layers is not None and "quantiles" in adata.layers:
                            row[f"{scorer_name}_quantile"] = float(adata.layers["quantiles"][adata_idx].sum())
            rows.append(row)

        df = pd.DataFrame(rows)
        os.makedirs("outputs", exist_ok=True)
        out_csv = "outputs/atlas_scn1a_summary.csv"
        df.to_csv(out_csv, index=False)
        print(f"Saved → {out_csv} ({len(df):,} rows)")

        print("\n=== Score distributions ===")
        for col in df.columns:
            if col.endswith("_score"):
                vals = df[col].dropna()
                if len(vals) > 0:
                    print(f"  {col}: mean={vals.mean():.4f}, std={vals.std():.4f}, "
                          f"top 1% threshold={vals.quantile(0.99):.4f}")

    print("\nDONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
