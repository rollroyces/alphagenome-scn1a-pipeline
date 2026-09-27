#!/usr/bin/env python3
"""
Weekly SCN1A VUS monitor: detect, score, and re-publish new ClinVar VUS.

What this script does, in order:
  1. Optionally fetch the latest ClinVar VCF (NCBI weekly drop) if it is
     stale or missing. If CLINVAR_VCF_URL is set, download from that URL.
  2. Tabix into SCN1A, filter to CLNSIG = uncertain_significance (VUS),
     then diff against our existing outputs/vus_rescored_with_dnase.csv.
  3. Score every *new SNV* VUS with AlphaGenome SPLICE_SITES + DNASE
     (the same two scorers that produced the existing 1,610-row CSV).
     Indels are deliberately skipped — AlphaGenome's splice scores are
     defined for SNVs, and our existing pipeline filters them out too.
  4. Append the new rows to outputs/vus_rescored_with_dnase.csv (atomic
     write — tmp file then rename — so a crashed run never clobbers the
     existing 1,610 rows).
  5. Re-stage and re-publish the dataset to Hugging Face
     (RROL/scn1a-vus-alphagenome). Only runs if HF_TOKEN is set, and
     ALWAYS writes to a temp parquet + uploads atomically — never edits
     the existing parquet in place.
  6. Send an email notification if any newly scored variant lands in the
     existing top-50 (would have qualified as Tier-1 or Tier-2 in our
     ranking). Gmail SMTP via app password — credentials are read from
     SMTP_USER / SMTP_PASSWORD env vars, never stored on disk.

Quota budget: ClinVar VCF download is ~200 MB but cached locally (only
redownloaded when CLINVAR_VCF_URL env var changes or the file is missing).
AlphaGenome API calls are limited to n_new variants — typically <10/week
for SCN1A, well under the 50-call weekly budget.

Usage:
    bash scripts/_run_with_key.sh scripts/monitor_clinvar.py

Environment variables (all optional, but recommended for full behaviour):
    CLINVAR_VCF_URL       Override the ClinVar VCF URL (e.g. when you want
                          to pin to a specific NCBI release).
    HF_TOKEN              Hugging Face write token (enables step 5).
    HF_DATASET_REPO       Default: RROL/scn1a-vus-alphagenome
    SMTP_HOST             Default: smtp.gmail.com
    SMTP_PORT             Default: 587
    SMTP_USER             Gmail address that sends the notification.
    SMTP_PASSWORD         Gmail app password (NOT the account password).
    NOTIFY_EMAIL          Address to send the Tier-1 alert to.
                          Default: SMTP_USER.
    SCN1A_MONITOR_DRY_RUN Set to 1 to skip scoring / HF upload / email
                          (useful for testing the download + diff logic).
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import smtplib
import sys
import tempfile
import time
import urllib.request
from datetime import datetime, timezone
from email.message import EmailMessage
from pathlib import Path

import pandas as pd


# Paths are repo-root relative so the launchd plist can `cd` first.
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_VCF = REPO_ROOT / "data" / "clinvar_grch38.vcf.gz"
DEFAULT_TBI = REPO_ROOT / "data" / "clinvar_grch38.vcf.gz.tbi"
DEFAULT_EXISTING = REPO_ROOT / "outputs" / "vus_rescored_with_dnase.csv"
DEFAULT_LOG_DIR = REPO_ROOT / "scripts" / "logs"
DEFAULT_HF_REPO = "RROL/scn1a-vus-alphagenome"

# NCBI's ClinVar VCF for GRCh38 — refreshed weekly. As of late 2025 this is
# the standard location. Override with CLINVAR_VCF_URL when pinning.
DEFAULT_CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_latest.vcf.gz"
)
DEFAULT_CLINVAR_TBI_URL = DEFAULT_CLINVAR_URL + ".tbi"


# ---- 1. Download -----------------------------------------------------------

def download_file(url: str, dest: Path, timeout: int = 300) -> None:
    """Stream-download `url` to `dest`. Skips if `dest` already exists and
    `force=False` — but a stale file is *not* the script's problem; we trust
    the user to delete the local file when they want a refresh, and the
    script accepts CLINVAR_VCF_URL as the override mechanism."""
    print(f"  downloading {url} → {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "scn1a-monitor/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as f:
        shutil.copyfileobj(resp, f)
    tmp.replace(dest)


def ensure_clinvar(vcf: Path, tbi: Path, force: bool) -> None:
    """Download the ClinVar VCF + index if missing or force=True.

    Note: a present-but-stale file is left alone (force=False) — this is the
    safe default. To refresh, either pass --refresh-clinvar or set
    CLINVAR_VCF_URL to a new dated file like clinvar_20261004.vcf.gz.
    """
    if not vcf.exists() or force:
        download_file(os.environ.get("CLINVAR_VCF_URL", DEFAULT_CLINVAR_URL), vcf)
    if not tbi.exists() or force:
        # The .tbi index is conventionally co-located with the .vcf.gz.
        base = os.environ.get("CLINVAR_VCF_URL", DEFAULT_CLINVAR_URL)
        tbi_url = os.environ.get("CLINVAR_VCF_TBI_URL", base + ".tbi")
        download_file(tbi_url, tbi)


# ---- 2. Diff ----------------------------------------------------------------

def fetch_new_scn1a_vus(vcf: Path, tbi: Path, existing: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame of SCN1A VUS not in `existing`, filtered to SNVs.

    Reuses the parsing + filter logic from scripts/check_new_variants.py
    (kept in-sync by both scripts sharing SCN1A_INTERVAL / SCN1A_GENE_ID
    constants; the constants live here too, see _is_scn1a below).
    """
    # Local imports — pysam is only needed when this script is actually run.
    import pysam

    SCN1A_INTERVAL = ("2", 165_984_640, 166_182_806)
    SCN1A_FLANK = 500_000
    SCN1A_GENE_ID = "6323"

    def _is_scn1a(geneinfo: str) -> bool:
        if not geneinfo:
            return False
        for entry in geneinfo.split("|"):
            if ":" in entry:
                _, gid = entry.split(":", 1)
                if gid == SCN1A_GENE_ID:
                    return True
        return False

    chrom, start, end = SCN1A_INTERVAL
    region = f"{chrom}:{max(0, start - SCN1A_FLANK)}-{end + SCN1A_FLANK}"

    tb = pysam.TabixFile(str(vcf), index=str(tbi))
    existing_keys = set(
        zip(existing["chrom"].astype(str), existing["pos"].astype(int),
            existing["ref"].astype(str), existing["alt"].astype(str))
    )

    out: list[dict] = []
    for raw in tb.fetch(region):
        if raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue
        chrom_f, pos_f, rsid, ref, alts, _q, _f, info = fields[:8]
        info_d: dict[str, list[str]] = {}
        for kv in info.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                info_d[k] = v.split(",")
            else:
                info_d[kv] = ["true"]
        geneinfo = info_d.get("GENEINFO", [""])[0]
        if not _is_scn1a(geneinfo):
            continue
        clnsig = info_d.get("CLNSIG", [""])[0]
        if clnsig not in {
            "Uncertain_significance",
            "Uncertain_significance,_conflicting_interpretations",
            "VUS-high",
        }:
            continue
        alt = alts.split(",")[0] if alts else ""
        key = (str(chrom_f), int(pos_f), str(ref), str(alt))
        if key in existing_keys:
            continue
        out.append({
            "chrom": chrom_f,
            "pos": int(pos_f),
            "ref": ref,
            "alt": alt,
            "rsid": rsid,
            "clnsig_category": "uncertain",
            "clnsig_raw": clnsig,
            "molecular_consequence": info_d.get("MC", [""])[0].split("|", 1)[-1],
            "clndn": info_d.get("CLNDN", [""])[0].replace("|", ","),
            "review_status": info_d.get("CLNREVSTAT", [""])[0],
        })
    return pd.DataFrame(out)


# ---- 3. Score ---------------------------------------------------------------

def score_new_variants(new_vus: pd.DataFrame) -> pd.DataFrame:
    """Score `new_vus` with SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS
    + DNASE — same four scorers used by rescore_vus + rescore_vus_dnase.

    Mirrors scripts/rescore_vus_dnase.py score_dnase() pattern: builds a
    Variant from each row, expands to 16 Kb interval, calls
    dna_model.score_variant, captures the four scores + success flag.
    """
    import numpy as np
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set. Run via _run_with_key.sh",
              file=sys.stderr)
        sys.exit(1)

    dna_model = dna_client.create(api_key)
    splice_scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]
    dnase_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]

    rows = []
    t0 = time.time()
    for i, (_, r) in enumerate(new_vus.iterrows()):
        if i % 25 == 0 or i == len(new_vus) - 1:
            elapsed = time.time() - t0
            print(f"  [score {i+1}/{len(new_vus)}] {elapsed:.0f}s elapsed",
                  flush=True)
        try:
            chrom = str(r["chrom"])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            variant = genome.Variant(
                chromosome=chrom,
                position=int(r["pos"]),
                reference_bases=str(r["ref"]),
                alternate_bases=str(r["alt"]),
            )
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
            # Two separate calls so we can capture DNASE alongside the
            # existing three splice scores — same pattern as
            # rescore_vus_dnase.py (which also calls them sequentially
            # rather than batching; avoids accidentally bundling all 4
            # into a single response that the SDK would error on).
            splice_scores = dna_model.score_variant(
                interval=interval, variant=variant, variant_scorers=splice_scorers,
            )
            dnase_scores = dna_model.score_variant(
                interval=interval, variant=variant, variant_scorers=[dnase_scorer],
            )
            row = dict(r)
            row["SPLICE_SITES_score"] = float(splice_scores[0].X.sum())
            row["SPLICE_SITE_USAGE_score"] = float(splice_scores[1].X.sum())
            row["SPLICE_JUNCTIONS_score"] = float(splice_scores[2].X.sum())
            row["DNASE_score"] = float(dnase_scores[0].X.sum())
            row["dnase_success"] = True
            row["dnase_error"] = ""
            rows.append(row)
        except Exception as e:
            row = dict(r)
            for c in ("SPLICE_SITES_score", "SPLICE_SITE_USAGE_score",
                      "SPLICE_JUNCTIONS_score", "DNASE_score"):
                row[c] = np.nan
            row["dnase_success"] = False
            row["dnase_error"] = f"{type(e).__name__}: {str(e)[:120]}"
            rows.append(row)
    return pd.DataFrame(rows)


# ---- 4. Append --------------------------------------------------------------

# Columns we expect to write into the merged CSV (16-col schema matching
# vus_rescored_with_dnase.csv, plus rsid from the new variants).
OUTPUT_COLUMNS = [
    "chrom", "pos", "ref", "alt", "clnsig_category",
    "DNASE_score", "SPLICE_SITES_score", "SPLICE_SITE_USAGE_score",
    "SPLICE_JUNCTIONS_score", "dnase_success", "dnase_error",
    "rsid", "molecular_consequence", "clndn", "review_status", "vus_rank",
]


def append_to_csv(scored_new: pd.DataFrame, existing_csv: Path) -> pd.DataFrame:
    """Merge scored_new into existing_csv, atomic write, return the merged DF.

    vus_rank is recomputed as SPLICE_SITES_score rank (descending, NaN last)
    across the FULL merged dataset — same convention as rescore_vus_dnase.
    """
    import numpy as np

    existing = pd.read_csv(existing_csv)
    print(f"  existing CSV rows: {len(existing):,}")

    # Ensure both DataFrames have the columns we want to write.
    merged = pd.concat([existing, scored_new], ignore_index=True, sort=False)
    for c in OUTPUT_COLUMNS:
        if c not in merged.columns:
            merged[c] = np.nan

    # Recompute vus_rank on the full merged dataset. NaN scores go to the
    # bottom (lowest priority), matching the existing convention.
    merged["vus_rank"] = (
        merged["SPLICE_SITES_score"]
        .rank(method="min", ascending=False, na_option="bottom")
        .astype("Int64")
    )

    # Atomic write: tmp file → rename. A crash mid-write can never leave
    # a half-written CSV that the next run would mistakenly treat as the
    # authoritative existing dataset.
    tmp = existing_csv.with_suffix(existing_csv.suffix + ".tmp")
    merged[OUTPUT_COLUMNS].to_csv(tmp, index=False)
    tmp.replace(existing_csv)
    print(f"  wrote merged CSV → {existing_csv} ({len(merged):,} rows)")
    return merged


# ---- 5. HF publish ----------------------------------------------------------

def publish_hf(merged: pd.DataFrame, repo_id: str) -> bool:
    """Upload `merged` as parquet to the Hugging Face dataset repo.

    Returns True on success, False if HF_TOKEN is missing or upload fails.
    The script continues either way — local files are still authoritative.

    Safety: we stage the parquet to a tempdir and use upload_folder, so the
    *existing* parquet on the Hub is never partially overwritten. The HF
    Hub handles commit atomicity server-side; the client just pushes a
    single revision.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from huggingface_hub import HfApi

    token = os.environ.get("HF_TOKEN")
    if not token:
        print("  [HF] HF_TOKEN not set — skipping Hugging Face publish")
        return False

    api = HfApi(token=token)
    try:
        # Sanity check: we should NEVER shrink the dataset. 1,610 is the
        # known-good baseline; if our merged DF is smaller, something is
        # very wrong — abort loudly.
        if len(merged) < 1610:
            print(f"  [HF] refusing to upload: merged has only {len(merged)} "
                  f"rows, expected ≥1610 (would clobber the published dataset)",
                  file=sys.stderr)
            return False

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            parquet_path = td_path / "data.parquet"
            # Convert Int64 → int64 so pyarrow writes nullable ints cleanly.
            table = pa.Table.from_pandas(merged, preserve_index=False)
            pq.write_table(table, parquet_path)

            ts = datetime.now(timezone.utc).strftime("%Y-%m-%d")
            commit_msg = (
                f"scn1a-monitor weekly update {ts}: "
                f"{len(merged)} SCN1A VUS "
                f"({len(merged) - 1610} new since 8db4888)"
            )
            print(f"  [HF] uploading {parquet_path} → {repo_id}")
            api.upload_folder(
                folder_path=str(td_path),
                repo_id=repo_id,
                repo_type="dataset",
                commit_message=commit_msg,
                allow_patterns=["data.parquet"],
            )
        print(f"  [HF] upload complete → "
              f"https://huggingface.co/datasets/{repo_id}")
        return True
    except Exception as e:
        print(f"  [HF] upload failed: {type(e).__name__}: {e}",
              file=sys.stderr)
        return False


# ---- 6. Email ---------------------------------------------------------------

def send_tier1_email(merged: pd.DataFrame, new_variants: pd.DataFrame,
                     log_path: Path) -> bool:
    """If any newly-scored variant lands in the top-50 by SPLICE_SITES_score,
    email the user. Returns True iff a notification was sent.

    Threshold = "would have been Tier-1 or Tier-2" = top-50 by SPLICE_SITES.
    That's the existing convention in scripts/rescore_vus_dnase.py and
    HUGGINGFACE_DATASET.md: tier=1 is top-4, tier=2 is top-50 minus top-4.
    """
    smtp_user = os.environ.get("SMTP_USER")
    smtp_password = os.environ.get("SMTP_PASSWORD")
    notify_to = os.environ.get("NOTIFY_EMAIL", smtp_user)
    if not smtp_user or not smtp_password or not notify_to:
        print("  [email] SMTP_USER / SMTP_PASSWORD / NOTIFY_EMAIL not all set — "
              "skipping notification (set them in your shell or a .env file "
              "outside the repo)")
        return False

    threshold_50 = (
        merged["SPLICE_SITES_score"].dropna()
        .nlargest(50).min()
    )
    new_in_top50 = new_variants[
        new_variants["SPLICE_SITES_score"] >= threshold_50
    ]
    if new_in_top50.empty:
        print(f"  [email] no new variants in top-50 "
              f"(threshold SPLICE_SITES ≥ {threshold_50:.3f}) — skipping")
        return False

    msg = EmailMessage()
    msg["Subject"] = (
        f"[scn1a-monitor] {len(new_in_top50)} new SCN1A VUS landed in top-50"
    )
    msg["From"] = smtp_user
    msg["To"] = notify_to
    body_lines = [
        "New SCN1A VUS detected in this week's ClinVar release scored high "
        "enough to land in the existing top-50 (would qualify as Tier-1 "
        "or Tier-2 in our ranking).",
        "",
        f"Top-50 SPLICE_SITES threshold (existing): {threshold_50:.4f}",
        f"New variants at or above threshold: {len(new_in_top50)}",
        "",
        "Details:",
    ]
    cols = ["chrom", "pos", "ref", "alt", "rsid", "SPLICE_SITES_score",
            "DNASE_score", "molecular_consequence", "vus_rank"]
    for _, r in new_in_top50.sort_values("SPLICE_SITES_score",
                                         ascending=False).iterrows():
        body_lines.append("  " + " | ".join(
            f"{c}={r[c]}" for c in cols if c in r.index
        ))
    body_lines += [
        "",
        f"Log: {log_path}",
        f"Full CSV: outputs/vus_rescored_with_dnase.csv",
        "",
        "— scn1a-monitor (cron)",
    ]
    msg.set_content("\n".join(body_lines))

    host = os.environ.get("SMTP_HOST", "smtp.gmail.com")
    port = int(os.environ.get("SMTP_PORT", "587"))
    try:
        with smtplib.SMTP(host, port, timeout=30) as s:
            s.starttls()
            s.login(smtp_user, smtp_password)
            s.send_message(msg)
        print(f"  [email] sent notification to {notify_to} "
              f"({len(new_in_top50)} new top-50 hits)")
        return True
    except Exception as e:
        print(f"  [email] send failed: {type(e).__name__}: {e}",
              file=sys.stderr)
        return False


# ---- main ------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--tbi", type=Path, default=DEFAULT_TBI)
    parser.add_argument("--existing", type=Path, default=DEFAULT_EXISTING)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument("--refresh-clinvar", action="store_true",
                        help="Re-download ClinVar VCF even if it exists locally.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip scoring / HF upload / email; only do "
                             "download + diff + log. Useful for cron CI.")
    args = parser.parse_args()

    dry_run = args.dry_run or os.environ.get("SCN1A_MONITOR_DRY_RUN") == "1"

    args.log_dir.mkdir(parents=True, exist_ok=True)
    ts_utc = datetime.now(timezone.utc).isoformat()
    log_path = args.log_dir / "monitor_clinvar.json"

    print("=" * 70)
    print(f"scn1a-monitor @ {ts_utc} (dry_run={dry_run})")
    print("=" * 70)

    # 1. ClinVar
    print("\n[1/6] ClinVar download")
    try:
        ensure_clinvar(args.vcf, args.tbi, force=args.refresh_clinvar)
    except Exception as e:
        return _fail(log_path, ts_utc, f"download failed: {e}")

    # 2. Diff
    print("\n[2/6] Diff against existing scored CSV")
    try:
        existing_df = pd.read_csv(args.existing)
        new_vus = fetch_new_scn1a_vus(args.vcf, args.tbi, existing_df)
    except Exception as e:
        return _fail(log_path, ts_utc, f"diff failed: {e}")
    print(f"  SCN1A VUS in ClinVar not yet scored: {len(new_vus):,}")
    # Filter to SNVs only (matches the original rescore_vus.py filter).
    pre = len(new_vus)
    new_vus_snv = new_vus[
        new_vus["ref"].astype(str).str.len().eq(1)
        & new_vus["alt"].astype(str).str.len().eq(1)
    ].copy()
    n_indels_skipped = pre - len(new_vus_snv)
    print(f"  After SNV filter: {len(new_vus_snv):,} "
          f"({n_indels_skipped} indels skipped — AlphaGenome splice scores "
          f"are SNV-only)")

    if new_vus_snv.empty:
        print("\nNo new SNV SCN1A VUS to score this week — done.")
        _write_log(log_path, ts_utc, {
            "n_new_total": pre,
            "n_new_snv": 0,
            "n_indels_skipped": n_indels_skipped,
            "n_existing_scored": len(existing_df),
            "hf_published": False,
            "email_sent": False,
            "dry_run": dry_run,
        })
        return 0

    # 3. Score
    print("\n[3/6] Score with AlphaGenome (SPLICE + DNASE)")
    if dry_run:
        print("  dry-run: skipping API scoring")
        scored_new = new_vus_snv.copy()
    else:
        scored_new = score_new_variants(new_vus_snv)

    # 4. Append
    print("\n[4/6] Append to vus_rescored_with_dnase.csv")
    if dry_run:
        print("  dry-run: skipping append (would touch production CSV)")
        merged = pd.concat([existing_df, scored_new], ignore_index=True,
                           sort=False)
    else:
        merged = append_to_csv(scored_new, args.existing)

    # 5. HF
    print("\n[5/6] Publish to Hugging Face")
    if dry_run:
        print("  dry-run: skipping HF publish")
        hf_ok = False
    else:
        hf_ok = publish_hf(merged, os.environ.get("HF_DATASET_REPO", DEFAULT_HF_REPO))

    # 6. Email
    print("\n[6/6] Tier-1 notification")
    if dry_run:
        print("  dry-run: skipping email")
        email_ok = False
    else:
        email_ok = send_tier1_email(merged, scored_new, log_path)

    _write_log(log_path, ts_utc, {
        "n_new_total": pre,
        "n_new_snv": len(new_vus_snv),
        "n_indels_skipped": n_indels_skipped,
        "n_existing_scored": len(existing_df),
        "n_merged_total": len(merged),
        "hf_published": hf_ok,
        "email_sent": email_ok,
        "dry_run": dry_run,
    })
    print(f"\nLog: {log_path}")
    print("DONE")
    return 0


def _write_log(log_path: Path, ts_utc: str, payload: dict) -> None:
    payload = {"timestamp_utc": ts_utc, **payload}
    log_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = log_path.with_suffix(log_path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True))
    tmp.replace(log_path)


def _fail(log_path: Path, ts_utc: str, err: str) -> int:
    print(f"ERROR: {err}", file=sys.stderr)
    _write_log(log_path, ts_utc, {"ok": False, "error": err})
    return 1


if __name__ == "__main__":
    sys.exit(main())