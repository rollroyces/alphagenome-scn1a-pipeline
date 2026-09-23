#!/usr/bin/env python3
"""
Experiment 011: ISM (in-silico mutagenesis) on the 4 SCN1A Tier-1 candidate
variants using the DNASE scorer.

Background
----------
The four SCN1A Tier-1 candidate variants (rs801806, rs4293437, rs801809,
rs2847163) were ranked at the top of outputs/vus_top_candidates.csv by
SPLICE_SITES_score. Three of them are annotated ClinVar
"splice_acceptor_variant" / "splice_donor_variant" consequences (the
fourth sits next to a known splice acceptor in the donor-acceptor pair);
they are prime candidates for the AlphaGenome splicing-driven Dravet story.

Experiment 004 already proved at n=30+30 per gene (DMD, CFTR, SCN1A) that
the fraction of DNASE ISM magnitude concentrated within +/-5bp of a
variant reliably distinguishes pathogenic from benign splice-region SNVs
(combined Mann-Whitney p=7.99e-7, r=+0.467). Here we ask the question:
"if we point AlphaGenome at the 4 Tier-1 candidates and ask it to ISM a
128 bp window centered on each one, what does it actually look like?"

The result is interpretable only in light of Exp 004's benchmark, since
n=4 has no statistical power on its own. The expected finding — if the
Tier-1 set is enriched for true splice-disrupting variants — is that the
+/-5bp concentration fraction is large (>0.5) and the maximum-effect
position falls very close to the variant.

Approach
--------
- 4 Tier-1 SCN1A variants from outputs/vus_top_candidates.csv
  (rs801806, rs4293437, rs801809, rs2847163).
- ISM window: 64 bp each side of the variant (128 bp total).
- Sequence length: 16,384 bp context (matches Exp 001-004).
- Scorer: ONLY DNASE CenterMaskScorer (the validated signal from Exp 004).
- Live API call to dna_model.score_ism_variants (NOT atlas.query_variants).
- Per-variant aggregation: average the 3 alt-allele rows, then read the
  central +/-5bp and +/-15bp concentration fractions from the resulting
  (128, 4) magnitude vector.
- Output: metrics.csv, 4 .npy matrices, 4-row heatmap figure, README.

Usage:
    bash scripts/_run_with_key.sh \
        research_notebook/experiments/011_ism_tier1/ism_dnase_tier1.py
"""

from __future__ import annotations

import csv
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# ---- Configuration --------------------------------------------------------

ISM_WINDOW = 64  # bp on each side of the variant (matches Exp 004)
CONTEXT = dna_client.SEQUENCE_LENGTH_16KB  # 16,384 bp context
MODALITY = "DNASE"

OUTPUT_DIR = Path("research_notebook/experiments/011_ism_tier1")
FIG_DIR = OUTPUT_DIR / "figures"

# Tier-1 SCN1A candidates from outputs/vus_top_candidates.csv
# (chrom, pos, ref, alt, rsid, molecular_consequence)
TIER1 = [
    ("chr2", 166041471, "T", "A", "rs801806", "splice_acceptor_variant"),
    ("chr2", 166073671, "C", "G", "rs4293437", "splice_acceptor_variant"),
    ("chr2", 166043700, "A", "G", "rs801809", "splice_donor_variant"),
    ("chr2", 166043701, "C", "A", "rs2847163", "splice_donor_variant"),
]

# Only DNASE this time — replicate the Exp 004 DNASE signature
SCORER = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]


# ---- Helpers (reuse / match Exp 004 conventions) --------------------------

def normalize_variant(v):
    """Return (chrom, pos, ref, alt) from a Variant object or 'chr:pos:ref>alt' string."""
    if isinstance(v, str):
        m = re.match(r"(.+):(\d+):([ACGT])>([ACGT])", v)
        if not m:
            raise ValueError(f"Cannot parse variant string: {v}")
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return str(v.chromosome), int(v.position), str(v.reference_bases), str(v.alternate_bases)


def run_ism(dna_model, chrom: str, pos: int, ref: str, alt: str, scorer):
    """Run ISM on one variant. Return list of (variant, total_score) per (position, alt)."""
    interval = genome.Interval(chrom, pos - CONTEXT // 2, pos + CONTEXT // 2)
    ism_interval = genome.Interval(chrom, pos - ISM_WINDOW, pos + ISM_WINDOW)

    variant_scores = dna_model.score_ism_variants(
        interval=interval,
        ism_interval=ism_interval,
        variant_scorers=[scorer],
    )

    results = []
    for variant_scores_list in variant_scores:
        for ann in variant_scores_list:
            # The API stores the Variant object under uns['variant'] (see
            # dna_client._construct_anndata_from_proto line ~356). Exp 004's
            # template uses uns['variant'] and it works.
            variant_obj = ann.uns.get("variant", "")
            total = float(np.asarray(ann.X).sum())
            results.append((variant_obj, total))
    return results


def build_ism_matrix(variant_results, center_pos):
    """Construct (n_positions, 4) ISM matrix from per-(position, alt) results.

    Ref base column stays zero (ISM only mutates to the other 3 bases).
    """
    bases = ["A", "C", "G", "T"]
    base_idx = {b: i for i, b in enumerate(bases)}
    n_positions = ISM_WINDOW * 2

    matrix = np.zeros((n_positions, 4), dtype=np.float32)
    for variant_obj, score in variant_results:
        _, pos, ref, alt = normalize_variant(variant_obj)
        rel_pos = pos - (center_pos - ISM_WINDOW)
        if 0 <= rel_pos < n_positions and alt in base_idx:
            matrix[rel_pos, base_idx[alt]] = score
    return matrix


def compute_concentration_metrics(ism_matrix):
    """Compute how concentrated ISM magnitude is near the variant (center)."""
    n_positions = ism_matrix.shape[0]
    center = n_positions // 2
    magnitude = np.abs(ism_matrix).sum(axis=1)  # sum over alt bases per position
    total = float(magnitude.sum())
    if total == 0:
        return {
            "total": 0.0,
            "fraction_within_5bp": 0.0,
            "fraction_within_15bp": 0.0,
            "concentration_5bp": 0.0,
            "max_position_relative_to_variant": -999,
            "max_value": 0.0,
        }
    within_5bp = float(magnitude[center - 5: center + 6].sum())
    within_15bp = float(magnitude[center - 15: center + 16].sum())
    max_pos = int(np.argmax(magnitude))
    max_val = float(magnitude[max_pos])
    return {
        "total": total,
        "fraction_within_5bp": within_5bp / total,
        "fraction_within_15bp": within_15bp / total,
        "concentration_5bp": within_5bp / (11 * magnitude.mean()) if magnitude.mean() > 0 else 0.0,
        "max_position_relative_to_variant": max_pos - center,
        "max_value": max_val,
    }


def classify_pattern(max_pos_rel: int, frac_5bp: float, consequence: str) -> str:
    """Heuristic splice-region classification from the ISM shape.

    - 'splice_donor_like'   : max-effect falls at +1..+6 bp (donor site is at the
                             5' splice junction, just downstream of the exon).
                             Strong +/-5bp concentration.
    - 'splice_acceptor_like': max-effect falls at -2..-1 bp (acceptor site is
                             at the 3' splice junction, just upstream of the
                             exon). Strong +/-5bp concentration.
    - 'broad_regulatory'    : max-effect far from variant, or weak +/-5bp
                             concentration: looks like a general regulatory
                             effect rather than a splice-site hit.
    """
    if frac_5bp < 0.30:
        return "broad_regulatory"
    if max_pos_rel <= -1:
        return "splice_acceptor_like"
    if max_pos_rel >= 1:
        return "splice_donor_like"
    return "splice_site_like_central"


# ---- Per-variant driver ---------------------------------------------------

def run_one(dna_model, chrom, pos, ref, alt, rsid, consequence):
    t0 = time.time()
    results = run_ism(dna_model, chrom, pos, ref, alt, SCORER)
    matrix = build_ism_matrix(results, pos)
    if matrix.shape != (ISM_WINDOW * 2, 4):
        raise ValueError(f"Got matrix shape {matrix.shape}, expected ({ISM_WINDOW * 2}, 4)")
    m = compute_concentration_metrics(matrix)
    classification = classify_pattern(m["max_position_relative_to_variant"],
                                      m["fraction_within_5bp"],
                                      consequence)
    m.update({
        "rsid": rsid,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "consequence": consequence,
        "classification": classification,
        "runtime_s": time.time() - t0,
        "success": True,
    })
    return m, matrix


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Creating dna_model client...")
    dna_model = dna_client.create(api_key)

    rows = []
    matrices = {}
    for chrom, pos, ref, alt, rsid, consequence in TIER1:
        print(f"\n--- {rsid}  {chrom}:{pos}:{ref}>{alt}  ({consequence}) ---")
        t0 = time.time()
        try:
            metrics, matrix = run_one(dna_model, chrom, pos, ref, alt, rsid, consequence)
            rows.append(metrics)
            matrices[rsid] = matrix
            np.save(OUTPUT_DIR / f"ism_dnase_{rsid}.npy", matrix)
            print(f"    +/-5bp frac={metrics['fraction_within_5bp']:.3f}, "
                  f"+/-15bp frac={metrics['fraction_within_15bp']:.3f}, "
                  f"max@{metrics['max_position_relative_to_variant']:+d}bp, "
                  f"total={metrics['total']:.1f} ({metrics['runtime_s']:.0f}s) -> "
                  f"{classification_label(metrics)}")
        except Exception as e:
            print(f"    FAILED: {type(e).__name__}: {str(e)[:200]}")
            rows.append({
                "rsid": rsid, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "consequence": consequence,
                "total": float("nan"), "fraction_within_5bp": float("nan"),
                "fraction_within_15bp": float("nan"), "concentration_5bp": float("nan"),
                "max_position_relative_to_variant": -999, "max_value": float("nan"),
                "classification": "FAILED",
                "runtime_s": time.time() - t0, "success": False,
                "error": f"{type(e).__name__}: {str(e)[:200]}",
            })

    # ---------- Write metrics.csv (4 rows, one per candidate) -----------------
    metric_cols = ["rsid", "chrom", "pos", "ref", "alt", "consequence",
                   "total", "fraction_within_5bp", "fraction_within_15bp",
                   "concentration_5bp", "max_position_relative_to_variant",
                   "max_value", "classification", "runtime_s"]
    with open(OUTPUT_DIR / "metrics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=metric_cols, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"\nSaved -> {OUTPUT_DIR / 'metrics.csv'} ({len(rows)} rows)")

    # ---------- Generate the 4-row heatmap ----------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rsids = [t[4] for t in TIER1]
        n_rows = len(rsids)
        n_pos = ISM_WINDOW * 2
        x = np.arange(-ISM_WINDOW, ISM_WINDOW)  # relative-to-variant positions

        # Per-position magnitude (sum across 3 alt bases), then log1p for
        # visualization. We add a tiny floor so positions with score ~0
        # render as deep blue rather than disappearing.
        fig, axes = plt.subplots(n_rows, 1, figsize=(11, 2.2 * n_rows),
                                  sharex=True)
        if n_rows == 1:
            axes = [axes]
        for ax, rsid in zip(axes, rsids):
            M = matrices.get(rsid)
            if M is None:
                ax.text(0.5, 0.5, f"{rsid}: no data",
                        ha="center", va="center", transform=ax.transAxes)
                continue
            mag = np.abs(M).sum(axis=1)
            log_mag = np.log1p(mag)
            im = ax.imshow(log_mag[np.newaxis, :], aspect="auto",
                           cmap="viridis",
                           extent=[-ISM_WINDOW, ISM_WINDOW, 0, 1])
            ax.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.9)
            # overlay markers at +/-5 and +/-15
            for w_, style in [(5, ":"), (15, ":")]:
                ax.axvline(-w_, color="white", linestyle=style, linewidth=0.6, alpha=0.4)
                ax.axvline(w_, color="white", linestyle=style, linewidth=0.6, alpha=0.4)
            ax.set_yticks([])
            ax.set_ylabel(rsid, rotation=0, ha="right", va="center", fontsize=10)
            # annotate +/-5bp concentration in the corner
            row = next(r for r in rows if r["rsid"] == rsid)
            ax.set_title(
                f"{rsid}  {row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}  "
                f"| frac+/-5bp={row['fraction_within_5bp']:.2f}  "
                f"frac+/-15bp={row['fraction_within_15bp']:.2f}  "
                f"max@{int(row['max_position_relative_to_variant']):+d}bp  "
                f"-> {row['classification']}",
                fontsize=8.5, loc="left")
        axes[-1].set_xlabel("Position relative to variant (bp)")
        fig.suptitle("ISM DNASE magnitude — SCN1A Tier-1 candidates",
                     fontsize=12, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        out = FIG_DIR / "tier1_ism_heatmap.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"Saved -> {out}")
    except Exception as e:
        print(f"Heatmap generation FAILED: {type(e).__name__}: {e}")

    print(f"\nDone. All outputs in {OUTPUT_DIR}/")
    return 0


def classification_label(m):
    return m["classification"]


if __name__ == "__main__":
    sys.exit(main())