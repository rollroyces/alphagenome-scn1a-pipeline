#!/usr/bin/env python3
"""
Query AlphaGenome Atlas for SCN1A — incremental CSV writer.

Smarter than atlas_scn1a_chunked.py: writes one row per variant to a CSV
as soon as each chunk completes. No AnnData accumulation in memory.
Survives crashes — you can re-run and it appends only new chunks.

Output: outputs/atlas_scn1a_summary.csv (one row per variant × allele)
Columns: chrom, pos, ref, alt, SPLICE_SITES_score, SPLICE_SITE_USAGE_score,
         SPLICE_JUNCTIONS_score, _quantile variants where available

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/atlas_scn1a_csv.py
"""

from __future__ import annotations

import csv
import os
import sys
import time

from alphagenome.atlas import atlas
from alphagenome.data import genome


SCN1A_CHROM = "chr2"
SCN1A_START = 165_984_640
SCN1A_END = 166_182_806
SCN1A_FLANK = 50_000
CHUNK_SIZE = 20_000  # smaller chunks = lower peak memory

SPLICING_SCORERS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
OUTPUT_CSV = "outputs/atlas_scn1a_summary.csv"


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    # Monkey-patch gRPC for large responses
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
    n_chunks = ((total_span - 1) // CHUNK_SIZE) + 1
    print(f"SCN1A region: {SCN1A_CHROM}:{start:,}-{end:,} ({total_span:,} bp)")
    print(f"Chunk size: {CHUNK_SIZE:,} bp → {n_chunks} chunks")
    print(f"Estimated variants: ~{total_span * 3:,}")

    os.makedirs("outputs", exist_ok=True)
    os.makedirs("data/atlas_progress", exist_ok=True)

    # Track which chunks are already done (for resumability)
    progress_dir = "data/atlas_progress"
    done_chunks = set()
    for fname in os.listdir(progress_dir):
        if fname.startswith("chunk_") and fname.endswith(".done"):
            try:
                idx = int(fname.split("_")[1].split(".")[0])
                done_chunks.add(idx)
            except ValueError:
                pass

    if done_chunks:
        print(f"Resuming — {len(done_chunks)} chunks already done")

    # Set up CSV
    file_exists = os.path.exists(OUTPUT_CSV)
    csv_file = open(OUTPUT_CSV, "a", newline="")
    fieldnames = [
        "chrom", "pos", "ref", "alt",
        "SPLICE_SITES_score", "SPLICE_SITE_USAGE_score", "SPLICE_JUNCTIONS_score",
        "SPLICE_SITES_quantile", "SPLICE_SITE_USAGE_quantile", "SPLICE_JUNCTIONS_quantile",
    ]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    if not file_exists:
        writer.writeheader()

    atlas_client = atlas.create(api_key)
    total_variants_written = 0
    t0 = time.time()

    for chunk_idx in range(1, n_chunks + 1):
        if chunk_idx in done_chunks:
            print(f"[{chunk_idx}/{n_chunks}] already done, skipping")
            continue

        chunk_start = start + (chunk_idx - 1) * CHUNK_SIZE
        chunk_end = min(chunk_start + CHUNK_SIZE, end)
        interval = genome.Interval(SCN1A_CHROM, chunk_start, chunk_end)
        print(f"\n[{chunk_idx}/{n_chunks}] {interval} ({interval.width:,} bp)")

        try:
            scores_by_scorer = atlas_client.query_interval(
                interval=interval,
                requested_scorers=SPLICING_SCORERS,
                max_workers=4,
                progress_bar=False,
            )

            # Pick the scorer with the most rows as the iteration index
            primary = max(scores_by_scorer.values(), key=lambda a: a.shape[0])
            n_variants = primary.shape[0]
            print(f"  Got {len(scores_by_scorer)} scorers, primary has {n_variants:,} variants")

            # Build per-variant rows
            chunk_rows = 0
            for i in range(n_variants):
                variant = primary.obs.iloc[i]["variant"]
                row = {
                    "chrom": variant.chromosome,
                    "pos": variant.position,
                    "ref": variant.reference_bases,
                    "alt": variant.alternate_bases,
                }
                for scorer_name, adata in scores_by_scorer.items():
                    if adata.obs is None or "variant" not in adata.obs.columns:
                        continue
                    # Find the row in this scorer matching (pos, ref, alt)
                    obs = adata.obs
                    mask = (
                        (obs["variant"].apply(lambda v: v.position) == variant.position)
                        & (obs["variant"].apply(lambda v: v.reference_bases) == variant.reference_bases)
                        & (obs["variant"].apply(lambda v: v.alternate_bases) == variant.alternate_bases)
                    )
                    if mask.any():
                        match_idx = mask.idxmax()
                        adata_idx = adata.obs.index.get_loc(match_idx)
                        row[f"{scorer_name}_score"] = float(adata.X[adata_idx].sum())
                        if adata.layers is not None and "quantiles" in adata.layers:
                            row[f"{scorer_name}_quantile"] = float(adata.layers["quantiles"][adata_idx].sum())
                writer.writerow(row)
                chunk_rows += 1

            csv_file.flush()
            os.sync = lambda: None  # noop
            with open(os.path.join(progress_dir, f"chunk_{chunk_idx:03d}.done"), "w") as f:
                f.write(f"{chunk_start}-{chunk_end}\n")
            total_variants_written += chunk_rows

            elapsed = time.time() - t0
            print(f"  ✓ Wrote {chunk_rows:,} variants to CSV")
            print(f"  Total: {total_variants_written:,} variants, elapsed {elapsed:.0f}s")

        except Exception as e:
            print(f"  ✗ Failed: {type(e).__name__}: {e}")
            print(f"  Recreating client and continuing...")
            try:
                del atlas_client
            except Exception:
                pass
            atlas_client = atlas.create(api_key)
            continue

    csv_file.close()
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print(f"DONE — {total_variants_written:,} variants written in {elapsed:.0f}s")
    print(f"Output: {OUTPUT_CSV}")
    print(f"To re-run safely: delete data/atlas_progress/*.done markers to redo chunks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
