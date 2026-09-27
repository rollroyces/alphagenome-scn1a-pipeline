#!/usr/bin/env python3
"""
Reproducibility smoke test for the SCN1A / AlphaGenome pipeline.

What it does
------------
1. Picks 5 SCN1A SNV variants at random (seeded, so the test is deterministic)
   from the existing `outputs/benchmark_scn1a_live_api_raw.csv` output of the
   full 591-variant benchmark run.
2. Re-scores each of the 5 variants against the live AlphaGenome API
   (`dna_model.score_variant` with SPLICE_SITES, SPLICE_SITE_USAGE,
   SPLICE_JUNCTIONS scorers, 16 Kb context).
3. Compares each re-scored value to the previously-recorded value in the raw
   CSV.  Variant scoring is deterministic for the same model + sequence +
   interval, so the two should match within a tight numeric tolerance.
4. Computes the 5-variant AUPRC and compares to the full-benchmark AUPRC
   (0.9833) — at n=5 we only require AUPRC > 0.5 (passes if all the picked
   variants are scored as before; the relative ranking must be preserved).
5. Prints PASS / FAIL with details, exits 0 on PASS, 1 on FAIL.

This burns 5 live API calls (3 scorers × 5 variants ≈ a few seconds). It is
intentionally cheap — the harness is for "did the install / API key / pipeline
still work", not for re-running the full benchmark.

Usage
-----
    export ALPHAGENOME_API_KEY=...
    python scripts/smoke_test.py

or, via the launcher that pulls the key from .alphagenome_key:

    bash scripts/_run_with_key.sh scripts/smoke_test.py

Environment variables
---------------------
SMOKE_TOLERANCE   : max allowed absolute difference per scorer
                    (default 1e-3 — same model = bit-identical in practice)
SMOKE_N_VARIANTS  : how many variants to re-score (default 5)
SMOKE_SEED        : RNG seed for picking the variants (default 42)
"""

from __future__ import annotations

import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_BENCH = REPO_ROOT / "outputs" / "benchmark_scn1a_live_api_raw.csv"
EXPECTED_RESULTS = REPO_ROOT / "outputs" / "benchmark_scn1a_results.csv"

SCORER_NAMES = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
N_VARIANTS = int(os.environ.get("SMOKE_N_VARIANTS", "5"))
SEED = int(os.environ.get("SMOKE_SEED", "42"))
TOLERANCE = float(os.environ.get("SMOKE_TOLERANCE", "1e-3"))

# Full-benchmark AUPRC for SPLICE_SITES (from benchmark_scn1a_results.csv).
# At n=5 we don't expect to *match* this — the point is the relative ranking
# of pathogenic vs benign is preserved.
FULL_AUPRC_SPLICE_SITES = 0.9833
# At n=5, a passing smoke test only requires the AUPRC to be > 0.5, i.e. the
# scoring is at least directionally correct.
MIN_AUPRC = 0.5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pick_variants(raw: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    """Pick n SNVs from the raw benchmark, biased to include both labels.

    We want the smoke test to be informative: pick roughly half pathogenic and
    half benign so the AUPRC > 0.5 check has a chance of triggering a real
    regression.
    """
    snvs = raw[
        (raw["ref"].str.len() == 1) & (raw["alt"].str.len() == 1)
        & raw["success"].astype(bool)
    ].copy()
    if len(snvs) < n:
        raise RuntimeError(
            f"Only {len(snvs)} successful SNVs in {RAW_BENCH}; need at least {n}."
        )

    pos = snvs[snvs["clnsig_category"] == "pathogenic"]
    neg = snvs[snvs["clnsig_category"] == "benign"]
    n_pos = min(n // 2 + 1, len(pos))
    n_neg = min(n - n_pos, len(neg))
    if n_pos + n_neg < n:
        # Fall back: top up from whichever side has surplus.
        shortfall = n - n_pos - n_neg
        if len(pos) - n_pos >= shortfall:
            n_pos += shortfall
        else:
            n_neg += shortfall

    rng = random.Random(seed)
    pos_sample = pos.sample(n=n_pos, random_state=rng.randrange(2**31))
    neg_sample = neg.sample(n=n_neg, random_state=rng.randrange(2**31))
    picked = pd.concat([pos_sample, neg_sample], ignore_index=True)
    return picked.reset_index(drop=True)


def _score_variants(dna_model, variants: pd.DataFrame) -> pd.DataFrame:
    """Re-score the chosen variants against the live API."""
    # Import inside the function so the module is importable even without the
    # alphagenome SDK (e.g. when building the Docker image's syntax-check).
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]

    out_rows = []
    for i, row in variants.iterrows():
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

        t0 = time.time()
        try:
            scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
        except Exception as e:  # pragma: no cover - depends on network
            print(f"  ! variant {i} ({row.get('rsid', '?')}) failed: {e}")
            out_rows.append({"index": i, "success": False, "error": str(e)})
            continue

        result = {"index": i, "success": True, "rsid": row.get("rsid", "")}
        for ann, name in zip(scores, SCORER_NAMES):
            result[f"{name}_score"] = float(ann.X.sum())
        result["elapsed_s"] = round(time.time() - t0, 2)
        out_rows.append(result)

    return pd.DataFrame(out_rows)


def _auprc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Tiny AUPRC for sanity check; defers to sklearn if available."""
    try:
        from sklearn.metrics import average_precision_score
        return float(average_precision_score(y_true, y_score))
    except Exception:
        # Fallback: trapezoidal PR curve (still correct, just slower).
        order = np.argsort(-y_score)
        yt = y_true[order]
        tp = np.cumsum(yt)
        fp = np.cumsum(1 - yt)
        recall = tp / tp[-1] if tp[-1] else np.zeros_like(tp, dtype=float)
        precision = tp / np.maximum(tp + fp, 1)
        # Trapezoidal integration under PR curve.
        ap = 0.0
        prev_r, prev_p = 0.0, 1.0
        for r, p in zip(recall, precision):
            ap += (r - prev_r) * (prev_p + p) / 2.0
            prev_r, prev_p = r, p
        return float(ap)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("FAIL: ALPHAGENOME_API_KEY environment variable is not set.")
        print("  Get a key at https://alphagenome.google/api and store it in")
        print("  .alphagenome_key, then run via scripts/_run_with_key.sh.")
        return 1

    if not RAW_BENCH.exists():
        print(f"FAIL: {RAW_BENCH} is missing. Cannot pick reference variants.")
        return 1

    raw = pd.read_csv(RAW_BENCH)
    print(f"[smoke] loaded {len(raw)} benchmark rows from {RAW_BENCH.name}")

    picked = _pick_variants(raw, N_VARIANTS, SEED)
    print(
        f"[smoke] picked {len(picked)} variants "
        f"({(picked['clnsig_category'] == 'pathogenic').sum()} pathogenic, "
        f"{(picked['clnsig_category'] == 'benign').sum()} benign) "
        f"with seed={SEED}"
    )

    # Import here so we fail fast on missing SDK without consuming API quota.
    from alphagenome.models import dna_client

    print("[smoke] creating dna_client …")
    dna_model = dna_client.create(api_key)

    print(f"[smoke] re-scoring {len(picked)} variants against the live API …")
    rescored = _score_variants(dna_model, picked)
    print(rescored.to_string(index=False))

    # ---- Compare against recorded values --------------------------------
    print(f"\n[smoke] tolerance per scorer: {TOLERANCE:g}")
    failures = []
    for _, scored in rescored.iterrows():
        if not scored["success"]:
            failures.append(f"variant {scored['index']} failed to score")
            continue
        original_row = picked.iloc[int(scored["index"])]
        for name in SCORER_NAMES:
            col = f"{name}_score"
            new = float(scored[col])
            old = float(original_row[col])
            diff = abs(new - old)
            status = "OK" if diff <= TOLERANCE else "FAIL"
            print(
                f"  {status}: variant {scored['index']} ({scored['rsid']}) "
                f"{name}: old={old:.6f} new={new:.6f} |Δ|={diff:.2e}"
            )
            if status == "FAIL":
                failures.append(
                    f"variant {scored['index']} {name} mismatch: "
                    f"|{new:.6f} - {old:.6f}|={diff:.2e} > {TOLERANCE:g}"
                )

    # ---- Mini-AUPRC sanity check -----------------------------------------
    labels = (picked["clnsig_category"] == "pathogenic").astype(int).to_numpy()
    scores = rescored["SPLICE_SITES_score"].to_numpy()
    if len(np.unique(labels)) > 1 and not np.isnan(scores).any():
        mini_auprc = _auprc(labels, scores)
    else:
        mini_auprc = float("nan")
    print(
        f"\n[smoke] mini-AUPRC (SPLICE_SITES, n={N_VARIANTS}): {mini_auprc:.4f} "
        f"(full benchmark AUPRC={FULL_AUPRC_SPLICE_SITES:.4f})"
    )
    if not np.isnan(mini_auprc) and mini_auprc < MIN_AUPRC:
        failures.append(
            f"mini-AUPRC {mini_auprc:.4f} < {MIN_AUPRC} — ranking broken"
        )

    # ---- Verdict --------------------------------------------------------
    if failures:
        print("\n[smoke] FAIL")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\n[smoke] PASS — pipeline reproduces existing outputs within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
