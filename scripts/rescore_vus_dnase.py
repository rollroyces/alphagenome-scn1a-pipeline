#!/usr/bin/env python3
"""
Experiment 009: re-score all 1,610 SCN1A VUS with DNASE in addition to
the existing splicing scores, then compute a joint ranking.

Inputs:
- outputs/vus_rescored.csv (1,610 VUS already scored on SPLICE_SITES,
  SPLICE_SITE_USAGE, SPLICE_JUNCTIONS — produced by scripts/rescore_vus.py)

Outputs (all in outputs/):
- vus_rescored_with_dnase.csv — 1,610 rows, all splicing scores + DNASE_score
- vus_top_by_dnase.csv        — top 30 by DNASE
- vus_top_combined.csv        — top 30 by mean(normalized_splicing, normalized_DNASE)
- vus_high_impact_with_dnase.csv — copy of vus_high_impact_with_gnomad.csv
  augmented with the new DNASE_score column (does NOT modify the original)

Usage:
    bash scripts/_run_with_key.sh python scripts/rescore_vus_dnase.py

Reuses the score_variant pattern from scripts/rescore_vus.py — same client,
same interval (16 Kb), same exception → NaN behavior. Throughput on Apple
Silicon is ~1.5 calls/sec → ~18 min for 1,610 variants.
"""

from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


VUS_INPUT = "outputs/vus_rescored.csv"
HIGH_IMPACT_INPUT = "outputs/vus_high_impact_with_gnomad.csv"

OUT_FULL = "outputs/vus_rescored_with_dnase.csv"
OUT_TOP_DNASE = "outputs/vus_top_by_dnase.csv"
OUT_TOP_COMBINED = "outputs/vus_top_combined.csv"
OUT_HIGH_DNASE = "outputs/vus_high_impact_with_dnase.csv"

TOP_N = 30


def score_dnase(vus: pd.DataFrame) -> pd.DataFrame:
    """Score each VUS with the DNASE scorer (live API)."""
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Run via _run_with_key.sh")
        sys.exit(1)

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    dnase_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]

    rows = []
    t0 = time.time()
    n_success = 0
    n_fail = 0

    for i, (_, row) in enumerate(vus.iterrows()):
        if i % 50 == 0 or i == len(vus) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(vus) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(vus)}] {elapsed:.0f}s elapsed, "
                  f"success={n_success}, fail={n_fail}, ETA {remaining:.0f}s",
                  flush=True)

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
                variant_scorers=[dnase_scorer],
            )
            dnase_score = float(scores[0].X.sum())
            result = dict(row)
            result["DNASE_score"] = dnase_score
            result["dnase_success"] = True
            result["dnase_error"] = ""
            n_success += 1
        except Exception as e:
            result = dict(row)
            result["DNASE_score"] = np.nan
            result["dnase_success"] = False
            result["dnase_error"] = f"{type(e).__name__}: {str(e)[:120]}"
            n_fail += 1

        rows.append(result)

    elapsed = time.time() - t0
    print(f"\nDNASE scoring done. {n_success}/{len(vus)} succeeded, "
          f"{n_fail} failed, total {elapsed:.0f}s")
    return pd.DataFrame(rows), elapsed, n_fail


def normalize_01(s: pd.Series) -> pd.Series:
    """Min-max normalize a series to [0, 1]. NaN-safe."""
    s = s.astype(float)
    valid = s.dropna()
    if len(valid) == 0:
        return s * np.nan
    lo, hi = valid.min(), valid.max()
    if hi == lo:
        return s * 0.0
    return (s - lo) / (hi - lo)


def main() -> int:
    print("=" * 70)
    print("Experiment 009 — DNASE re-scoring of SCN1A VUS")
    print("=" * 70)

    # 1. Load the existing VUS splice-scored file
    if not os.path.exists(VUS_INPUT):
        print(f"ERROR: {VUS_INPUT} not found. Run scripts/rescore_vus.py first.")
        return 1
    vus = pd.read_csv(VUS_INPUT)
    print(f"Loaded {len(vus):,} VUS from {VUS_INPUT}")
    required = {"chrom", "pos", "ref", "alt", "SPLICE_SITES_score"}
    missing = required - set(vus.columns)
    if missing:
        print(f"ERROR: missing required columns in VUS input: {missing}")
        return 1

    # 2. Score DNASE
    scored, runtime_s, n_fail = score_dnase(vus)
    n_scored = (~scored["DNASE_score"].isna()).sum()
    print(f"\nSummary: scored {n_scored:,}/{len(vus):,}, "
          f"failed {n_fail:,} ({n_fail / len(vus) * 100:.1f}%), "
          f"runtime {runtime_s:.0f}s")

    # Ensure clnsig_category column exists (the original vus_rescored.csv stores it
    # but the column may already have the categorical value 'uncertain_significance')
    if "clnsig_category" not in scored.columns:
        scored["clnsig_category"] = "uncertain_significance"

    # 3. Save the raw DNASE-augmented file
    out_cols = ["chrom", "pos", "ref", "alt", "clnsig_category", "DNASE_score",
                "SPLICE_SITES_score", "SPLICE_SITE_USAGE_score",
                "SPLICE_JUNCTIONS_score", "dnase_success", "dnase_error",
                "rsid", "molecular_consequence", "clndn", "review_status",
                "vus_rank"]
    for c in out_cols:
        if c not in scored.columns:
            scored[c] = np.nan
    scored_sorted = scored.sort_values("vus_rank", na_position="last")
    scored_sorted[out_cols].to_csv(OUT_FULL, index=False)
    print(f"Saved full DNASE-augmented file → {OUT_FULL} "
          f"({len(scored_sorted):,} rows)")

    # 4. Top 30 by DNASE
    dnase_valid = scored.dropna(subset=["DNASE_score"]).copy()
    dnase_valid = dnase_valid.sort_values("DNASE_score", ascending=False)
    dnase_top = dnase_valid.head(TOP_N).copy()
    dnase_top["dnase_rank"] = range(1, len(dnase_top) + 1)
    cols_top = ["dnase_rank", "chrom", "pos", "ref", "alt", "rsid",
                "DNASE_score", "SPLICE_SITES_score", "molecular_consequence",
                "clndn", "review_status"]
    for c in cols_top:
        if c not in dnase_top.columns:
            dnase_top[c] = ""
    dnase_top[cols_top].to_csv(OUT_TOP_DNASE, index=False)
    print(f"Saved top {TOP_N} by DNASE → {OUT_TOP_DNASE}")

    # 5. Joint ranking: mean(normalized SPLICE_SITES, normalized DNASE).
    # SPLICE_SITES was the best single discriminator in rescore_vus.py.
    joint = scored.dropna(subset=["DNASE_score", "SPLICE_SITES_score"]).copy()
    joint["splice_norm"] = normalize_01(joint["SPLICE_SITES_score"])
    joint["dnase_norm"] = normalize_01(joint["DNASE_score"])
    joint["combined_score"] = (joint["splice_norm"] + joint["dnase_norm"]) / 2.0
    joint = joint.sort_values("combined_score", ascending=False)
    joint_top = joint.head(TOP_N).copy()
    joint_top["combined_rank"] = range(1, len(joint_top) + 1)

    # 6. Splicing-only top-30 by SPSEN rankings on the SAME filtered set
    #    (so the comparison is apples-to-apples)
    spl_only = scored.dropna(subset=["SPLICE_SITES_score"]).copy()
    spl_only_top = spl_only.sort_values(
        "SPLICE_SITES_score", ascending=False
    ).head(TOP_N)
    spl_only_ids = set(
        zip(spl_only_top["chrom"], spl_only_top["pos"],
            spl_only_top["ref"], spl_only_top["alt"])
    )
    joint_top_ids = set(
        zip(joint_top["chrom"], joint_top["pos"],
            joint_top["ref"], joint_top["alt"])
    )
    overlap = spl_only_ids & joint_top_ids
    only_in_splicing = spl_only_ids - joint_top_ids
    only_in_combined = joint_top_ids - spl_only_ids

    cols_joint = ["combined_rank", "chrom", "pos", "ref", "alt", "rsid",
                  "combined_score", "splice_norm", "dnase_norm",
                  "SPLICE_SITES_score", "DNASE_score",
                  "molecular_consequence", "clndn", "review_status"]
    for c in cols_joint:
        if c not in joint_top.columns:
            joint_top[c] = ""
    joint_top[cols_joint].to_csv(OUT_TOP_COMBINED, index=False)
    print(f"Saved top {TOP_N} by combined → {OUT_TOP_COMBINED}")

    print(f"\n=== Ranking comparison (top {TOP_N}) ===")
    print(f"  Splicing-only top-30: {len(spl_only_ids)} variants")
    print(f"  Combined top-30:      {len(joint_top_ids)} variants")
    print(f"  Overlap:              {len(overlap)}")
    print(f"  Only in splicing:     {len(only_in_splicing)}")
    print(f"  Only in combined:     {len(only_in_combined)}")
    if len(overlap) >= 20:
        verdict = "REDUNDANT — DNASE adds little new signal vs splicing alone."
    elif len(overlap) <= 10:
        verdict = "COMPLEMENTARY — DNASE pulls in many variants splicing missed."
    else:
        verdict = "MIXED — partial overlap; DNASE reshuffles but isn't fully redundant."
    print(f"  Verdict: {verdict}")

    print(f"\n=== Top 10 by combined ranking ===")
    display = joint_top.head(10)[
        ["combined_rank", "rsid", "chrom", "pos", "ref", "alt",
         "SPLICE_SITES_score", "DNASE_score", "combined_score",
         "molecular_consequence"]
    ]
    print(display.to_string(index=False))

    # 7. Augment high-impact with gnomAD with DNASE (NEW file, original untouched)
    if os.path.exists(HIGH_IMPACT_INPUT):
        hi = pd.read_csv(HIGH_IMPACT_INPUT)
        key_cols = ["chrom", "pos", "ref", "alt"]
        dnase_map = scored[["chrom", "pos", "ref", "alt", "DNASE_score",
                            "dnase_success"]].rename(
            columns={"DNASE_score": "DNASE_score_new"}
        )
        merged = hi.merge(dnase_map, on=key_cols, how="left")
        # Use DNASE_score if the original column already exists, else fall back
        if "DNASE_score" in merged.columns:
            merged["DNASE_score"] = merged["DNASE_score"].fillna(
                merged["DNASE_score_new"]
            )
        else:
            merged["DNASE_score"] = merged["DNASE_score_new"]
        merged = merged.drop(columns=["DNASE_score_new"])
        merged.to_csv(OUT_HIGH_DNASE, index=False)
        n_with_dnase = merged["DNASE_score"].notna().sum()
        print(f"\nSaved high-impact + DNASE → {OUT_HIGH_DNASE} "
              f"({n_with_dnase}/{len(merged)} have DNASE_score)")
    else:
        print(f"\nWARNING: {HIGH_IMPACT_INPUT} not found; "
              f"skipped {OUT_HIGH_DNASE}")

    # Stash a comparison summary next to the run
    summary_path = "outputs/vus_combined_vs_splicing.json"
    import json
    summary = {
        "n_scored": int(n_scored),
        "n_failed": int(n_fail),
        "n_total": int(len(vus)),
        "runtime_seconds": round(runtime_s, 1),
        "overlap_top30": int(len(overlap)),
        "only_in_splicing_top30": int(len(only_in_splicing)),
        "only_in_combined_top30": int(len(only_in_combined)),
        "verdict": verdict,
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved comparison summary → {summary_path}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Outputs:")
    print(f"  {OUT_FULL}            ← all VUS with DNASE + splicing scores")
    print(f"  {OUT_TOP_DNASE}       ← top {TOP_N} by DNASE only")
    print(f"  {OUT_TOP_COMBINED}    ← top {TOP_N} by combined (key artifact)")
    print(f"  {OUT_HIGH_DNASE}      ← 61 high-impact candidates + DNASE")
    return 0


if __name__ == "__main__":
    sys.exit(main())