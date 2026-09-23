#!/usr/bin/env python3
"""
Experiment 014: Re-score the 4 SCN1A Tier-1 candidates using the RNA_SEQ scorer
to predict brain expression effect for Carvill lab prioritization.

Background
----------
Carvill lab can test 1-2 SCN1A variants per year in their minigene splicing
assay. We have 4 Tier-1 candidates (splice_donor_variant or splice_acceptor_variant
in ClinVar) and need concrete data to pick the most impactful 1-2.

Splicing assays measure RNA-level outcome, and SCN1A is highly brain-expressed.
RNA-seq brain tissue signals are the most relevant proxy for predicted haploinsufficiency
severity (more brain RNA-seq change -> more likely the allele reduces functional
transcript in the tissue where SCN1A actually matters).

Inputs
------
- outputs/vus_high_impact_candidates.csv (4 splice donor/acceptor variants)

Outputs
-------
- outputs/tier1_rnaseq_features.csv       (4 rows, per-track scores + aggregates)
- research_notebook/experiments/014_tier1_rnaseq/README.md
- research_notebook/experiments/014_tier1_rnaseq/tier1_rnaseq_per_track.csv (long
  format: one row per (variant, track) for downstream plotting)

Usage
-----
    bash scripts/_run_with_key.sh python scripts/rescore_tier1_rnaseq.py

NOTE: We use `variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']` directly
because this SDK version does not expose `dna_client.Scorer(preset_name=...)`.
The result is identical: an RNA_SEQ GeneMaskLFCScorer returning 371 RNA-seq
tracks. Same call pattern as scripts/rescore_vus_dnase.py.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# --- Inputs / outputs --------------------------------------------------------

# The Tier-1 splice donor/acceptor variants live in vus_high_impact_candidates.csv
# (the README/top_candidates.csv file holds a different splice ranking).
HIGH_IMPACT_INPUT = "outputs/vus_high_impact_candidates.csv"

TIER1_RSID_SET = {"801806", "4293437", "801809", "2847163"}

OUT_FEATURES = "outputs/tier1_rnaseq_features.csv"
OUT_PER_TRACK = "research_notebook/experiments/014_tier1_rnaseq/tier1_rnaseq_per_track.csv"

# Track-level "brain" filter: GTEx brain tissues + ENCODE brain biosamples.
# AlphaGenome's GTEx RNA-seq tracks label brain regions with the 'Brain_'
# prefix (e.g. 'Brain_Cerebellum', 'Brain_Frontal_Cortex_BA9'). Anything
# starting with 'Brain_' is treated as brain.

# Substring tokens to catch brain biosamples that don't go through GTEx (e.g.
# ENCODE brain biosamples from brain donors, or biosample_name == 'brain' / 'Brain').
# Tokens are matched in lowercase. Order matters: specific phrases must come
# before single words that could match elsewhere (e.g. "cortex" alone matches
# "kidney cortex" which is NOT brain).
BRAIN_BIOSAMPLE_TOKENS = (
    "frontal cortex",
    "prefrontal cortex",
    "anterior cingulate cortex",
    "occipital lobe",
    "parietal lobe",
    "temporal lobe",
    "dorsolateral prefrontal",
    "motor neuron",
    "hippocamp",
    "thalam",
    "amygdala",
    "ganglia",
    "cerebellum",
    "cerebellar",
    "neural progenitor",
    "neuronal stem",
    "neurosphere",
    "glutamatergic",
    "gabaergic",
    "dopaminergic neuron",
    "astrocyte",
    "microglia",
    "oligodendrocyte",
    "brain",
)


def is_brain_track(row: pd.Series) -> bool:
    """Decide whether a track row represents a brain tissue.

    Combines two signals:
      1. gtex_tissue starts with 'Brain_' (GTEx naming convention used by
         AlphaGenome RNA-seq tracks: 'Brain_Cerebellum', 'Brain_Cortex', etc.).
      2. biosample_name contains a brain-related token (catches ENCODE brain
         biosamples and SCN1A-relevant cell types such as iPSC-derived neurons).
    """
    gtex = str(row.get("gtex_tissue", "") or "")
    if gtex.startswith("Brain_"):
        return True
    bio = str(row.get("biosample_name", "") or "").lower()
    if any(tok in bio for tok in BRAIN_BIOSAMPLE_TOKENS):
        return True
    return False


def score_one_variant(dna_model, chrom, pos, ref, alt) -> anndata.AnnData:
    """Run dna_model.score_variant for a single variant with RNA_SEQ scorer."""
    if not str(chrom).startswith("chr"):
        chrom = "chr" + str(chrom)
    variant = genome.Variant(
        chromosome=chrom,
        position=int(pos),
        reference_bases=str(ref),
        alternate_bases=str(alt),
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
    rnase_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
    scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=[rnase_scorer],
    )
    return scores[0]


def aggregate_per_track(adata, rsid: str) -> pd.DataFrame:
    """Return a long-format per-track DataFrame for one variant's AnnData."""
    var_df = adata.var.copy()
    var_df["rsid"] = rsid

    # Extract the SCN1A row from X.
    scn1a_mask = adata.obs["gene_name"].astype(str).str.upper().eq("SCN1A")
    if not scn1a_mask.any():
        # Fallback: pick any gene; SCN1A should be the closest hit but log clearly.
        sys.stderr.write(
            f"[WARN] rs{rsid}: no SCN1A row in obs (genes: "
            f"{adata.obs['gene_name'].astype(str).str.upper().tolist()[:5]}...). "
            "Using row 0 as a placeholder.\n"
        )
        row_idx = 0
    else:
        row_idx = int(np.argmax(scn1a_mask.values))
    var_df["RNA_SEQ_score"] = np.asarray(adata.X[row_idx, :]).astype(float)
    var_df["is_brain"] = var_df.apply(is_brain_track, axis=1)
    return var_df


def main() -> int:
    print("=" * 70)
    print("Experiment 014 — RNA_SEQ re-scoring of 4 SCN1A Tier-1 candidates")
    print("=" * 70)

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Run via _run_with_key.sh")
        return 1
    if not os.path.exists(HIGH_IMPACT_INPUT):
        print(f"ERROR: {HIGH_IMPACT_INPUT} not found.")
        return 1

    candidates = pd.read_csv(HIGH_IMPACT_INPUT)
    # Tier-1 = splice_donor_variant or splice_acceptor_variant
    tier1 = candidates[
        candidates["molecular_consequence"].isin(
            ["splice_donor_variant", "splice_acceptor_variant"]
        )
    ].copy()
    tier1 = tier1[tier1["rsid"].astype(str).isin(TIER1_RSID_SET)].copy()
    if len(tier1) != 4:
        print(
            f"ERROR: expected 4 Tier-1 variants, got {len(tier1)}. "
            f"Found: {tier1['rsid'].astype(str).tolist()}"
        )
        return 1
    # Stable order: rs801806, rs4293437, rs801809, rs2847163
    order = {r: i for i, r in enumerate(
        ["801806", "4293437", "801809", "2847163"]
    )}
    tier1 = tier1.sort_values(
        by="rsid", key=lambda s: s.astype(str).map(order)
    ).reset_index(drop=True)
    print(f"Loaded {len(tier1)} Tier-1 candidates from {HIGH_IMPACT_INPUT}:")
    print(
        tier1[["rsid", "chrom", "pos", "ref", "alt", "molecular_consequence"]]
        .to_string(index=False)
    )

    print("\nCreating dna_model client...")
    dna_model = dna_client.create(api_key)

    per_track_long = []
    summary_rows = []
    t0 = time.time()

    for i, row in tier1.iterrows():
        rsid = str(row["rsid"])
        t_run = time.time()
        print(
            f"\n[{i+1}/{len(tier1)}] rs{rsid} "
            f"chr{int(row['chrom'])}:{int(row['pos'])} "
            f"{row['ref']}>{row['alt']} ({row['molecular_consequence']})"
        )
        try:
            adata = score_one_variant(
                dna_model,
                chrom=row["chrom"],
                pos=int(row["pos"]),
                ref=row["ref"],
                alt=row["alt"],
            )
            elapsed = time.time() - t_run
            n_tracks = int(adata.X.shape[1])
            per_track_df = aggregate_per_track(adata, rsid)
            per_track_long.append(per_track_df)

            total_score = float(per_track_df["RNA_SEQ_score"].sum())
            abs_score = per_track_df["RNA_SEQ_score"].abs()
            max_abs = float(abs_score.max())
            n_sig = int((abs_score > 1.0).sum())
            brain = per_track_df[per_track_df["is_brain"]]
            brain_total = float(brain["RNA_SEQ_score"].sum())
            brain_max_abs = float(brain["RNA_SEQ_score"].abs().max())
            n_brain = int(len(brain))
            n_brain_sig = int((brain["RNA_SEQ_score"].abs() > 1.0).sum())

            print(
                f"  ok ({elapsed:.1f}s). "
                f"total={total_score:+.3f} "
                f"brain_total={brain_total:+.3f} "
                f"max|track|={max_abs:.3f} "
                f"n_sig={n_sig}/{n_tracks} "
                f"brain_n={n_brain} (sig {n_brain_sig})"
            )

            summary_rows.append({
                "rsid": rsid,
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "molecular_consequence": row["molecular_consequence"],
                "RNA_SEQ_total": total_score,
                "RNA_SEQ_brain_total": brain_total,
                "RNA_SEQ_max_abs": max_abs,
                "RNA_SEQ_n_significant_tracks": n_sig,
                "RNA_SEQ_n_total_tracks": n_tracks,
                "RNA_SEQ_n_brain_tracks": n_brain,
                "RNA_SEQ_n_brain_significant": n_brain_sig,
                "RNA_SEQ_brain_max_abs": brain_max_abs,
            })
        except Exception as e:
            elapsed = time.time() - t_run
            print(f"  FAILED ({elapsed:.1f}s): {type(e).__name__}: {e}")
            summary_rows.append({
                "rsid": rsid,
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "molecular_consequence": row["molecular_consequence"],
                "RNA_SEQ_total": np.nan,
                "RNA_SEQ_brain_total": np.nan,
                "RNA_SEQ_max_abs": np.nan,
                "RNA_SEQ_n_significant_tracks": 0,
                "RNA_SEQ_n_total_tracks": 0,
                "RNA_SEQ_n_brain_tracks": 0,
                "RNA_SEQ_n_brain_significant": 0,
                "RNA_SEQ_brain_max_abs": np.nan,
            })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_FEATURES, index=False)
    print(f"\nSaved per-variant features -> {OUT_FEATURES}")

    long_df = pd.concat(per_track_long, ignore_index=True)
    long_df = long_df[
        [
            "rsid",
            "name",
            "strand",
            "biosample_name",
            "biosample_type",
            "gtex_tissue",
            "Assay title",
            "is_brain",
            "RNA_SEQ_score",
        ]
    ].sort_values(["rsid", "RNA_SEQ_score"], ascending=[True, False])
    os.makedirs(os.path.dirname(OUT_PER_TRACK), exist_ok=True)
    long_df.to_csv(OUT_PER_TRACK, index=False)
    print(f"Saved long-form per-track scores -> {OUT_PER_TRACK} "
          f"({len(long_df)} rows)")

    elapsed_total = time.time() - t0
    print(f"\nTotal runtime: {elapsed_total:.0f}s")

    # Print ranking preview to stdout (used in summary back to parent).
    print("\n=== Ranking by |brain_total| (higher = stronger predicted brain effect) ===")
    if summary["RNA_SEQ_brain_total"].notna().any():
        summary_ranked = summary.copy()
        summary_ranked["abs_brain_total"] = (
            summary_ranked["RNA_SEQ_brain_total"].abs()
        )
        summary_ranked = summary_ranked.sort_values(
            "abs_brain_total", ascending=False
        )
        print(
            summary_ranked[
                [
                    "rsid",
                    "RNA_SEQ_total",
                    "RNA_SEQ_brain_total",
                    "RNA_SEQ_brain_max_abs",
                    "RNA_SEQ_n_brain_significant",
                    "RNA_SEQ_n_brain_tracks",
                    "molecular_consequence",
                ]
            ].to_string(index=False)
        )
    print("\n=== DONE ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())