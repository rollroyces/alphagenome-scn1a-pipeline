#!/usr/bin/env python3
"""
Score a pre-extracted ClinVar benchmark for a single gene using AlphaGenome.

Reads a TSV from extract_clinvar_for_gene.py and scores each variant.
Writes per-gene raw scores CSV.

Usage:
    bash scripts/_run_with_key.sh scripts/score_gene_benchmark.py \
        --input outputs/clinvar_scn1a_benchmark.tsv \
        --output outputs/cross_disease_scn1a_raw.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


SCORER_NAMES = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]


def score_one(dna_model, chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Score a single variant. Returns dict with success flag and scores."""
    try:
        chrom_str = str(chrom)
        if not chrom_str.startswith("chr"):
            chrom_str = "chr" + chrom_str
        variant = genome.Variant(
            chromosome=chrom_str,
            position=int(pos),
            reference_bases=str(ref),
            alternate_bases=str(alt),
        )
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
        scorers = [
            variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in SCORER_NAMES
        ]
        scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
        )
        out = {"score_success": True}
        for ann, name in zip(scores, SCORER_NAMES):
            out[f"{name}_score"] = float(ann.X.sum())
        return out
    except Exception as e:
        return {
            "score_success": False,
            "score_error": f"{type(e).__name__}: {str(e)[:80]}",
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--gene", default=None, help="Override gene name in output")
    args = parser.parse_args()

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Use scripts/_run_with_key.sh")
        return 1

    df = pd.read_csv(args.input, sep="\t")
    gene = args.gene or df["gene"].iloc[0] if "gene" in df.columns else "unknown"
    print(f"Scoring {len(df)} variants for {gene}")
    print(f"  Positives: {(df['label'] == 'positive').sum()}")
    print(f"  Negatives: {(df['label'] == 'negative').sum()}")

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    rows = []
    t0 = time.time()
    n_success = 0
    n_fail = 0

    for i, row in df.iterrows():
        if i % 25 == 0 or i == len(df) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(df) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(df)}] {elapsed:.0f}s elapsed, "
                  f"success={n_success}, fail={n_fail}, ETA {remaining:.0f}s", flush=True)

        result = score_one(dna_model, row["chrom"], row["pos"], row["ref"], row["alt"])
        result.update({
            "gene": gene,
            "chrom": row["chrom"], "pos": row["pos"],
            "ref": row["ref"], "alt": row["alt"],
            "rsid": row.get("rsid", ""),
            "label": row["label"],
            "clnsig_category": row.get("clnsig_category", ""),
            "molecular_consequence": row.get("molecular_consequence", ""),
        })
        rows.append(result)
        if result["score_success"]:
            n_success += 1
        else:
            n_fail += 1

    elapsed = time.time() - t0
    print(f"\nDone. {n_success}/{len(df)} succeeded, {n_fail} failed, total {elapsed:.0f}s")

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out_df.to_csv(args.output, index=False)
    print(f"Saved → {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
