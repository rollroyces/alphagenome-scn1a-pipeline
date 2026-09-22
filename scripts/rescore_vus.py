#!/usr/bin/env python3
"""
Re-score ClinVar VUS (variants of uncertain significance) in SCN1A using
AlphaGenome's live API. Output: ranked candidate list of variants most likely
to be pathogenic based on splicing impact score.

This is the deliverable that goes to the Carvill / Sparber / Helbig labs.

Input:
- outputs/clinvar_scn1a.tsv (already exists, 5,276 SCN1A variants)

Output:
- outputs/vus_rescored.csv — every VUS with AlphaGenome score + rank
- outputs/vus_top_candidates.csv — top 50 candidates with annotations

Runtime: ~10–15 minutes for ~1,700 VUS

Usage:
    bash scripts/_run_with_key.sh  # exports key from .alphagenome_key
    python scripts/rescore_vus.py
"""

from __future__ import annotations

import csv
import os
import sys
import time

import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


SCN1A_OUTPUT = "outputs/clinvar_scn1a.tsv"
VUS_OUTPUT = "outputs/vus_rescored.csv"
TOP_OUTPUT = "outputs/vus_top_candidates.csv"


def load_vus() -> pd.DataFrame:
    """Load ClinVar SCN1A variants, filter to VUS."""
    df = pd.read_csv(SCN1A_OUTPUT, sep="\t")
    print(f"Loaded {len(df):,} total SCN1A variants from {SCN1A_OUTPUT}")

    # VUS = uncertain_significance category
    vus = df[df["clnsig_category"] == "uncertain"].copy()
    print(f"VUS (uncertain significance): {len(vus):,}")

    # Exclude variants that are non-SNVs (AlphaGenome splicing scores are
    # designed for SNVs and small indels). Skip indels and complex variants.
    # ClinVar doesn't have a clean filter; we'll skip multi-base alts.
    pre = len(vus)
    vus = vus[vus["ref"].str.len() == 1].copy()
    vus = vus[vus["alt"].str.len() == 1].copy()
    print(f"After filter to SNVs only: {len(vus):,} (dropped {pre - len(vus):,} indels)")

    return vus


def score_vus(vus: pd.DataFrame) -> pd.DataFrame:
    """Score each VUS with SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS."""
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Run via _run_with_key.sh")
        sys.exit(1)

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    # Use all three splicing scorers
    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]
    scorer_names = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]

    rows = []
    t0 = time.time()
    n_success = 0
    n_fail = 0

    for i, (_, row) in enumerate(vus.iterrows()):
        if i % 25 == 0 or i == len(vus) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(vus) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(vus)}] {elapsed:.0f}s elapsed, "
                  f"success={n_success}, fail={n_fail}, ETA {remaining:.0f}s")

        try:
            chrom = str(row["chrom"])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            variant = genome.Variant(
                chromosome=chrom,
                position=int(row["pos"]),
                reference_bases=str(row["ref"]),
                alternate_bases=str(row["alt"]),
            )
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
            scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
            result = dict(row)
            for ann, name in zip(scores, scorer_names):
                result[f"{name}_score"] = float(ann.X.sum())
            result["score_success"] = True
            n_success += 1
        except Exception as e:
            result = dict(row)
            result["score_success"] = False
            result["score_error"] = f"{type(e).__name__}: {str(e)[:100]}"
            n_fail += 1

        rows.append(result)

    elapsed = time.time() - t0
    print(f"\nDone. {n_success}/{len(vus)} succeeded, {n_fail} failed, total {elapsed:.0f}s")
    return pd.DataFrame(rows)


def rank_and_save(df: pd.DataFrame) -> None:
    """Rank variants by SPLICE_SITES score and save outputs."""
    # Filter to successful scores
    successful = df[df["score_success"] == True].copy()
    failed = df[df["score_success"] != True]
    print(f"\nSuccessful: {len(successful):,}, Failed: {len(failed):,}")

    if len(failed) > 0:
        print("Sample failures:")
        print(failed[["chrom", "pos", "ref", "alt", "score_error"]].head(5).to_string())

    # Rank by SPLICE_SITES score (best discriminator from benchmark)
    successful = successful.sort_values("SPLICE_SITES_score", ascending=False)
    successful["vus_rank"] = range(1, len(successful) + 1)

    # Save full ranked list
    os.makedirs("outputs", exist_ok=True)
    successful.to_csv(VUS_OUTPUT, index=False)
    print(f"\nSaved ranked VUS → {VUS_OUTPUT}")

    # Top 50 candidates with key fields
    top_cols = [
        "vus_rank", "chrom", "pos", "ref", "alt", "rsid",
        "SPLICE_SITES_score", "SPLICE_SITE_USAGE_score", "SPLICE_JUNCTIONS_score",
        "molecular_consequence", "clndn", "review_status",
    ]
    # Add columns if they exist
    for col in top_cols:
        if col not in successful.columns:
            print(f"  (missing column: {col})")

    available_cols = [c for c in top_cols if c in successful.columns]
    top50 = successful.head(50)[available_cols].copy()
    top50.to_csv(TOP_OUTPUT, index=False)
    print(f"Saved top 50 → {TOP_OUTPUT}")

    # Headline stats
    print("\n=== Top 10 VUS candidates (predicted splice impact) ===")
    display_cols = ["vus_rank", "chrom", "pos", "ref", "alt", "SPLICE_SITES_score", "molecular_consequence"]
    display_cols = [c for c in display_cols if c in successful.columns]
    print(successful.head(10)[display_cols].to_string(index=False))

    print("\n=== Score distribution of VUS ===")
    print(successful["SPLICE_SITES_score"].describe().to_string())

    # Compare with the 0.5 threshold from benchmark (high-impact = pathogenic-like)
    high_impact = successful[successful["SPLICE_SITES_score"] >= 0.5]
    print(f"\n★ High-impact VUS (score >= 0.5, pathogenic-like): {len(high_impact)}")
    print(f"  These are the priority candidates for experimental validation")

    # Save high-impact subset as separate file
    if len(high_impact) > 0:
        high_path = "outputs/vus_high_impact_candidates.csv"
        high_impact[available_cols].to_csv(high_path, index=False)
        print(f"  Saved → {high_path}")


def main() -> int:
    print("=" * 70)
    print("SCN1A VUS re-scoring with AlphaGenome")
    print("=" * 70)

    vus = load_vus()
    if len(vus) == 0:
        print("No VUS to score — done")
        return 0

    scored = score_vus(vus)
    rank_and_save(scored)

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Outputs:")
    print(f"  {VUS_OUTPUT}  ← all VUS with scores, ranked")
    print(f"  {TOP_OUTPUT}  ← top 50 candidates for outreach")
    print(f"  outputs/vus_high_impact_candidates.csv  ← score >= 0.5")
    print(f"\nNext:")
    print(f"  1. Review outputs/vus_high_impact_candidates.csv")
    print(f"  2. Cross-reference with literature for known candidates")
    print(f"  3. Prepare outreach email to Carvill / Sparber labs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
