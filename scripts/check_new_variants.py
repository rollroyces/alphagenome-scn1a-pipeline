#!/usr/bin/env python3
"""
Lightweight SCN1A VUS monitor: download the latest ClinVar VCF, count any
*new* SCN1A VUS that are NOT in our existing 1,610, and log the result.

This script makes ZERO AlphaGenome API calls — its job is to detect "is there
something new to score?" so a separate (heavier) weekly job can act on the
signal. Daily frequency is fine because all it does is parse a tabix-indexed
VCF (already cached locally after first run) and a 1,610-row CSV.

Usage:
    # default paths — repo-root relative
    python scripts/check_new_variants.py

    # explicit override (useful for cron / tests)
    python scripts/check_new_variants.py --vcf data/clinvar_grch38.vcf.gz \\
        --existing outputs/vus_rescored_with_dnase.csv \\
        --log scripts/logs/check_new_variants.json

Outputs:
    - Human-readable summary on stdout (n_new, n_total, n_existing).
    - JSON status file (--log) so downstream code / cron can diff between runs.

Exit codes:
    0  — ran cleanly (whether new variants were found or not)
    1  — VCF / index / existing-CSV missing or malformed
    2  — unexpected exception during parsing
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import pysam


# SCN1A gene interval (hg38) — kept in sync with scripts/clinvar_scn1a_local.py
SCN1A_INTERVAL = ("2", 165_984_640, 166_182_806)
SCN1A_FLANK = 500_000
SCN1A_GENE_ID = "6323"

UNCERTAIN_LABELS = {
    "Uncertain_significance",
    "Uncertain_significance,_conflicting_interpretations",
    "VUS-high",
}


def _is_scn1a(geneinfo: str) -> bool:
    """Check whether SCN1A (Gene ID 6323) appears in ClinVar's GENEINFO field.

    ClinVar separates multiple genes with a pipe, not a comma:
        GENEINFO=SCN1A:6323|LOC102724058:102724058
    """
    if not geneinfo:
        return False
    for entry in geneinfo.split("|"):
        if ":" in entry:
            _, gene_id = entry.split(":", 1)
            if gene_id == SCN1A_GENE_ID:
                return True
    return False


def _parse_record(line: str) -> dict | None:
    """Parse a single tabix-returned VCF line into a structured dict."""
    if line.startswith("#"):
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return None
    chrom, pos, rsid, ref, alts, _qual, _filt, info = fields[:8]
    info_dict: dict[str, list[str]] = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info_dict[k] = v.split(",")
        else:
            info_dict[kv] = ["true"]
    return {
        "chrom": chrom,
        "pos": int(pos),
        "rsid": rsid,
        "ref": ref,
        "alts": alts.split(","),
        "info": info_dict,
    }


def fetch_scn1a_vus(vcf_path: str, tbi_path: str) -> list[dict]:
    """Return every ClinVar SCN1A variant with CLNSIG = uncertain."""
    chrom, start, end = SCN1A_INTERVAL
    region = f"{chrom}:{max(0, start - SCN1A_FLANK)}-{end + SCN1A_FLANK}"
    tb = pysam.TabixFile(vcf_path, index=tbi_path)
    out: list[dict] = []
    for raw in tb.fetch(region):
        rec = _parse_record(raw)
        if rec is None:
            continue
        geneinfo = rec["info"].get("GENEINFO", [""])[0]
        if not _is_scn1a(geneinfo):
            continue
        clnsig = rec["info"].get("CLNSIG", [""])[0]
        if clnsig not in UNCERTAIN_LABELS:
            continue
        out.append({
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": rec["alts"][0] if rec["alts"] else "",
            "rsid": rec["rsid"],
            "clnsig": clnsig,
            "molecular_consequence": rec["info"].get("MC", [""])[0].split("|", 1)[-1],
        })
    return out


def load_existing_keys(existing_csv: str) -> set[tuple[str, int, str, str]]:
    """Return the set of (chrom, pos, ref, alt) tuples we have already scored.

    rsIDs are unreliable — ClinVar uses dbSNP rsIDs that change with dbSNP
    builds, but chrom/pos/ref/alt is stable for a given assembly. That tuple
    is what we use to detect "new" variants.
    """
    import pandas as pd

    df = pd.read_csv(existing_csv)
    return set(zip(df["chrom"].astype(str), df["pos"].astype(int),
                   df["ref"].astype(str), df["alt"].astype(str)))


def write_log(log_path: Path, payload: dict) -> None:
    """Write a JSON status file. Writes atomically (tmp → rename) so a crashed
    run never leaves a half-written file that the next monitor run would
    misinterpret as 'success'."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = log_path.with_suffix(log_path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True))
    tmp.replace(log_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vcf",
        default="data/clinvar_grch38.vcf.gz",
        help="Path to tabix-indexed ClinVar VCF (.vcf.gz). Default: %(default)s",
    )
    parser.add_argument(
        "--tbi",
        default="data/clinvar_grch38.vcf.gz.tbi",
        help="Path to VCF index (.vcf.gz.tbi). Default: %(default)s",
    )
    parser.add_argument(
        "--existing",
        default="outputs/vus_rescored_with_dnase.csv",
        help="Path to the existing scored CSV. New SCN1A VUS are detected by "
             "absent (chrom, pos, ref, alt) tuples. Default: %(default)s",
    )
    parser.add_argument(
        "--log",
        default="scripts/logs/check_new_variants.json",
        help="Where to write the JSON status file. Default: %(default)s",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the per-run summary on stdout.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.vcf) or not os.path.exists(args.tbi):
        print(f"ERROR: {args.vcf} or {args.tbi} not found.", file=sys.stderr)
        return 1
    if not os.path.exists(args.existing):
        print(f"ERROR: existing scored CSV {args.existing} not found.",
              file=sys.stderr)
        return 1

    ts = datetime.now(timezone.utc).isoformat()
    try:
        all_scn1a_vus = fetch_scn1a_vus(args.vcf, args.tbi)
        existing = load_existing_keys(args.existing)
    except Exception as e:
        err_payload = {
            "timestamp_utc": ts,
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
        }
        try:
            write_log(Path(args.log), err_payload)
        except Exception:
            pass  # best-effort
        print(f"ERROR: failed to parse ClinVar: {e}", file=sys.stderr)
        return 2

    new_keys: list[tuple[str, int, str, str]] = []
    for v in all_scn1a_vus:
        key = (str(v["chrom"]), int(v["pos"]), str(v["ref"]), str(v["alt"]))
        if key not in existing:
            new_keys.append(key)

    payload = {
        "timestamp_utc": ts,
        "ok": True,
        "vcf_path": args.vcf,
        "existing_csv": args.existing,
        "n_scn1a_vus_total": len(all_scn1a_vus),
        "n_existing_scored": len(existing),
        "n_new": len(new_keys),
        "new_variants": [
            {"chrom": c, "pos": p, "ref": r, "alt": a}
            for (c, p, r, a) in new_keys
        ],
    }
    write_log(Path(args.log), payload)

    if not args.quiet:
        print(f"[check_new_variants] {ts}")
        print(f"  ClinVar SCN1A VUS total: {len(all_scn1a_vus):,}")
        print(f"  Already scored:          {len(existing):,}")
        print(f"  NEW (not yet scored):    {len(new_keys):,}")
        for c, p, r, a in new_keys[:10]:
            print(f"    {c}:{p} {r}>{a}")
        if len(new_keys) > 10:
            print(f"    ... and {len(new_keys) - 10} more (see {args.log})")
    return 0


if __name__ == "__main__":
    sys.exit(main())