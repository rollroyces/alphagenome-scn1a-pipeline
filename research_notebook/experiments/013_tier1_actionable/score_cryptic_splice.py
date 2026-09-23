#!/usr/bin/env python3
"""
Experiment 013 — Fresh AlphaGenome score_variant on the 4 Tier-1 SCN1A
candidates with the explicit splicing scorers (SPLICE_SITES,
SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) at 16 KB context. The values in
outputs/vus_top_candidates.csv were originally taken from the cached
Atlas query; this re-runs the live model so the lab sees a focused,
self-contained cryptic-splice assessment for each candidate.

Usage:
    bash scripts/_run_with_key.sh \
        research_notebook/experiments/013_tier1_actionable/score_cryptic_splice.py

Output:
    research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.json
    research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.csv
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


# ---- Configuration --------------------------------------------------------

CONTEXT = dna_client.SEQUENCE_LENGTH_16KB  # 16,384 bp context
OUTPUT_DIR = Path("research_notebook/experiments/013_tier1_actionable")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Tier-1 SCN1A candidates
TIER1 = [
    ("chr2", 166041471, "T", "A", "rs801806", "splice_acceptor_variant"),
    ("chr2", 166073671, "C", "G", "rs4293437", "splice_acceptor_variant"),
    ("chr2", 166043700, "A", "G", "rs801809", "splice_donor_variant"),
    ("chr2", 166043701, "C", "A", "rs2847163", "splice_donor_variant"),
]

# Explicit splicing scorers (Atlas does not contain SPLICE_JUNCTIONS_ACTIVE)
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

    # Patch gRPC for large messages (same as Exp 011)
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

    print("Creating dna_model client...")
    model = dna_client.create(api_key)

    # Resolve scorer objects
    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    scorer_objs = []
    for s in SPLICING_SCORERS:
        if s not in recommended:
            print(f"⚠ Scorer {s} not in RECOMMENDED_VARIANT_SCORERS — skipping")
            continue
        scorer_objs.append(recommended[s])
        print(f"  Using scorer: {s}")
    print(f"Will score with {len(scorer_objs)} scorers\n")

    # Resolve variant objects
    variants = []
    for chrom, pos, ref, alt, rsid, mc in TIER1:
        v = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref,
                           alternate_bases=alt, name=rsid)
        variants.append((v, rsid, mc))

    # Score each variant on its own. Use variant.reference_interval.resize
    # to centre the 16KB window on the variant (matches Exp 011 pattern).
    results = []
    for v, rsid, mc in variants:
        t0 = time.time()
        print(f"Scoring {rsid} ({v.chromosome}:{v.position} {v.reference_bases}>{v.alternate_bases})...")
        try:
            interval = v.reference_interval.resize(CONTEXT)
            adatas = model.score_variant(
                interval=interval,
                variant=v,
                variant_scorers=scorer_objs,
                organism=dna_client.Organism.HOMO_SAPIENS,
            )
        except Exception as e:
            print(f"  ✗ Failed: {e}")
            results.append({
                'rsid': rsid, 'chrom': v.chromosome, 'pos': v.position,
                'ref': v.reference_bases, 'alt': v.alternate_bases,
                'molecular_consequence': mc, 'error': str(e),
            })
            continue
        dt = time.time() - t0
        # AnnData structure: adata.X is (1 variant x n_tracks). For SPLICE_SITES
        # and SPLICE_SITE_USAGE, .X holds alt - ref signed deltas; we report
        # the abs sum (cryptic splice disruption magnitude, the convention
        # used elsewhere in the project, e.g. outputs/vus_top_candidates.csv).
        # For SPLICE_JUNCTIONS, .X is the per-junction score for the alt allele
        # (a non-signed scorer).
        row = {
            'rsid': rsid, 'chrom': v.chromosome, 'pos': v.position,
            'ref': v.reference_bases, 'alt': v.alternate_bases,
            'molecular_consequence': mc, 'latency_s': round(dt, 1),
            'sequence_length': CONTEXT,
        }
        # Parse the requested_output name from the scorer's name string:
        # 'GeneMaskSplicingScorer(requested_output=SPLICE_SITES, width=None)' → 'SPLICE_SITES'
        # 'SpliceJunctionScorer()' → 'SPLICE_JUNCTIONS'
        import re
        def _short(s):
            m = re.search(r'requested_output=([A-Z_]+)', s)
            if m:
                return m.group(1)
            # SpliceJunctionScorer
            if 'SpliceJunction' in s:
                return 'SPLICE_JUNCTIONS'
            return s
        for adata, scorer_obj in zip(adatas, scorer_objs):
            scorer_name = _short(scorer_obj.name)
            x = np.asarray(adata.X)
            if x.ndim == 2:
                vals = x[0]
            else:
                vals = x
            row[f"{scorer_name}_sum"] = float(np.nansum(vals))
            row[f"{scorer_name}_abs_sum"] = float(np.nansum(np.abs(vals)))
            row[f"{scorer_name}_max_abs"] = float(np.nanmax(np.abs(vals)))
            row[f"{scorer_name}_n_tracks"] = int(vals.size)
        results.append(row)
        print(f"  ✓ {dt:.1f}s | SPLICE_SITES abs_sum={row.get('SPLICE_SITES_abs_sum'):.3f} "
              f"SPLICE_SITE_USAGE abs_sum={row.get('SPLICE_SITE_USAGE_abs_sum'):.3f} "
              f"SPLICE_JUNCTIONS abs_sum={row.get('SPLICE_JUNCTIONS_abs_sum'):.3f}")

    # Persist
    out_json = OUTPUT_DIR / "cryptic_splice_scores.json"
    out_csv = OUTPUT_DIR / "cryptic_splice_scores.csv"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)

    if results:
        fieldnames = list(results[0].keys())
        with open(out_csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in results:
                w.writerow(r)

    print(f"\nSaved {out_json}")
    print(f"Saved {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
