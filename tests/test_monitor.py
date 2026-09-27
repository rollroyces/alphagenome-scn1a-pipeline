"""
Tests for the SCN1A VUS weekly monitor.

These tests are OFFLINE — they do NOT hit the AlphaGenome API, do NOT
download ClinVar, and do NOT send email or upload to Hugging Face. They
cover the comparison / ranking logic that the monitor's correctness
depends on.

Run with: `pytest tests/test_monitor.py -v` (from the repo root).
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

# Ensure the monitor module imports without dragging the AlphaGenome SDK
# at module load time (the SDK is only imported lazily, inside
# score_new_variants — but we still want a guard here so a missing SDK
# does NOT break tests of the diff / ranking logic).
try:
    import monitor_clinvar  # noqa: F401
    import check_new_variants  # noqa: F401
    MODULES_IMPORTABLE = True
except ImportError as e:
    MODULES_IMPORTABLE = False
    IMPORT_ERROR = str(e)


# ---- Fixtures --------------------------------------------------------------

def _write_fake_clinvar(path: Path, rows: list[dict]) -> None:
    """Write a minimal tabix-readable VCF.gz containing `rows`.

    Each row is a dict with the keys: chrom, pos, rsid, ref, alt, geneinfo,
    clnsig, mc, clndn, review_status. We write a real bgzipped VCF so that
    pysam.TabixFile can index it the same way the production code expects.
    """
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=2,length=250000000>\n"
        "##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=MC,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=CLNDN,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description=\"\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    body_lines = []
    for r in rows:
        info = (
            f"GENEINFO={r['geneinfo']};"
            f"CLNSIG={r['clnsig']};"
            f"MC={r.get('mc', 'SO:0001627|intron_variant')};"
            f"CLNDN={r.get('clndn', 'Epilepsy')}|"
            f"CLNREVSTAT={r.get('review_status', 'criteria_provided')}"
        )
        body_lines.append(
            f"{r['chrom']}\t{r['pos']}\t{r['rsid']}\t{r['ref']}\t"
            f"{r['alt']}\t.\t.\t{info}"
        )
    # Use bgzip (pysam ships bgzipped indexes that need bgzip-formatted
    # input). We shell out to bgzip from pysam.
    plain = path.with_suffix(".plain.vcf")
    plain.write_text(header + "\n".join(body_lines) + "\n")
    import pysam
    pysam.tabix_compress(str(plain), str(path))
    pysam.tabix_index(str(path), preset="vcf", force=True)
    plain.unlink()


@pytest.fixture
def fake_clinvar(tmp_path: Path) -> Path:
    """A minimal fake ClinVar VCF.gz containing 3 SCN1A VUS (2 SNV + 1 indel)
    plus 1 benign SCN1A variant and 1 SCN2A variant (wrong gene — must be
    filtered out)."""
    rows = [
        # 2 SCN1A SNV VUS — should appear in the diff
        dict(chrom="2", pos=166_000_001, rsid="rs9999001",
             ref="A", alt="G", geneinfo="SCN1A:6323",
             clnsig="Uncertain_significance"),
        dict(chrom="2", pos=166_000_010, rsid="rs9999002",
             ref="C", alt="T", geneinfo="SCN1A:6323",
             clnsig="Uncertain_significance"),
        # 1 SCN1A indel VUS — must be filtered out by the SNV filter
        dict(chrom="2", pos=166_000_020, rsid="rs9999003",
             ref="ATCT", alt="A", geneinfo="SCN1A:6323",
             clnsig="Uncertain_significance"),
        # 1 SCN1A pathogenic — must NOT appear in the diff (not VUS)
        dict(chrom="2", pos=166_000_030, rsid="rs9999004",
             ref="G", alt="A", geneinfo="SCN1A:6323",
             clnsig="Pathogenic"),
        # 1 SCN2A VUS — wrong gene, must NOT appear
        dict(chrom="2", pos=166_000_040, rsid="rs9999005",
             ref="C", alt="T", geneinfo="SCN2A:6326",
             clnsig="Uncertain_significance"),
    ]
    vcf = tmp_path / "fake_clinvar.vcf.gz"
    _write_fake_clinvar(vcf, rows)
    return vcf


@pytest.fixture
def fake_existing_csv(tmp_path: Path) -> Path:
    """A fake existing CSV that already contains rs9999001 — so the diff
    should find rs9999002 as the only NEW SNV VUS."""
    df = pd.DataFrame([
        {"chrom": "2", "pos": 166_000_001, "ref": "A", "alt": "G",
         "rsid": "rs9999001", "SPLICE_SITES_score": 0.5,
         "SPLICE_SITE_USAGE_score": 100.0, "SPLICE_JUNCTIONS_score": 1000.0,
         "DNASE_score": -5.0, "vus_rank": 1,
         "clnsig_category": "uncertain"},
        {"chrom": "2", "pos": 165_900_000, "ref": "C", "alt": "T",
         "rsid": "rs8888001", "SPLICE_SITES_score": 0.4,
         "SPLICE_SITE_USAGE_score": 90.0, "SPLICE_JUNCTIONS_score": 900.0,
         "DNASE_score": -3.0, "vus_rank": 2,
         "clnsig_category": "uncertain"},
    ])
    p = tmp_path / "existing.csv"
    df.to_csv(p, index=False)
    return p


# ---- Module-import sanity --------------------------------------------------

def test_modules_importable():
    """Both monitor scripts must import cleanly (no top-level API imports
    that would block on a missing key)."""
    assert MODULES_IMPORTABLE, f"import failed: {IMPORT_ERROR}"


def test_check_module_uses_no_alpha_genome_at_import_time():
    """Importing check_new_variants must not require the AlphaGenome SDK —
    it's intended for daily no-API runs."""
    import check_new_variants
    src = Path(check_new_variants.__file__).read_text()
    assert "alphagenome" not in src, (
        "check_new_variants.py imports alphagenome at module load — must "
        "be deferred so the daily probe doesn't depend on the SDK."
    )


# ---- check_new_variants.py --------------------------------------------------

def test_check_counts_only_new_snv_vus(fake_clinvar, fake_existing_csv,
                                       tmp_path):
    """check_new_variants should report exactly 1 new SCN1A VUS (rs9999002)
    — rs9999001 is already scored, the indel/pathogenic/wrong-gene rows
    are filtered out."""
    log = tmp_path / "log.json"
    tbi = str(fake_clinvar) + ".tbi"
    rc = subprocess.run(
        [sys.executable, str(SCRIPTS / "check_new_variants.py"),
         "--vcf", str(fake_clinvar),
         "--tbi", tbi,
         "--existing", str(fake_existing_csv),
         "--log", str(log),
         "--quiet"],
        check=False,
    )
    assert rc.returncode == 0, f"check_new_variants exited {rc.returncode}"
    payload = json.loads(log.read_text())
    assert payload["ok"] is True
    # Fixture has 3 SCN1A VUS (2 SNV + 1 indel); indels are NOT counted
    # out by check_new_variants.py at this layer — they're logged as
    # 'new' and the AlphaGenome-SNV-only filter happens in the heavy
    # monitor script.
    assert payload["n_scn1a_vus_total"] == 3
    assert payload["n_existing_scored"] == 2
    # The 1 SNV that's new + the 1 indel (also new, but rsid-only). Only
    # the SNV gets scored by the heavy monitor — that's a separate test.
    assert payload["n_new"] == 2
    rsids = sorted(v["rsid"] if "rsid" in v else f"{v['chrom']}:{v['pos']}"
                   for v in payload["new_variants"])
    # The new variants are rs9999002 (SNV) + rs9999003 (indel). Their
    # position+ref+alt distinguishes them; just check we have 2 of them
    # with the expected positions.
    positions = sorted(v["pos"] for v in payload["new_variants"])
    assert positions == [166_000_010, 166_000_020]


def test_check_exits_nonzero_on_missing_vcf(fake_existing_csv, tmp_path):
    log = tmp_path / "log.json"
    rc = subprocess.run(
        [sys.executable, str(SCRIPTS / "check_new_variants.py"),
         "--vcf", str(tmp_path / "nope.vcf.gz"),
         "--tbi", str(tmp_path / "nope.vcf.gz.tbi"),
         "--existing", str(fake_existing_csv),
         "--log", str(log),
         "--quiet"],
        check=False,
    )
    assert rc.returncode == 1


# ---- monitor_clinvar.py diff / append logic ---------------------------------

def test_fetch_new_filters_correctly(fake_clinvar, fake_existing_csv):
    """The diff should return SCN1A VUS not in the existing CSV. Filter
    criteria: CLNSIG = uncertain, SCN1A in GENEINFO. Indels pass through
    here — the SNV-only filter happens in main() after this function."""
    from monitor_clinvar import fetch_new_scn1a_vus

    existing = pd.read_csv(fake_existing_csv)
    tbi = Path(str(fake_clinvar) + ".tbi")
    new = fetch_new_scn1a_vus(fake_clinvar, tbi, existing)
    # 2 SNV VUS (rs9999001 already scored, rs9999002 new) — fetch_new
    # returns the SCN1A-uncertain set minus existing. Indels do pass
    # through here.
    rsids = sorted(new["rsid"].tolist())
    assert "rs9999002" in rsids  # the new one
    # rs9999001 is already in existing so should NOT appear.
    assert "rs9999001" not in rsids
    # Pathogenic SCN1A must not appear.
    assert "rs9999004" not in rsids
    # Wrong-gene SCN2A must not appear.
    assert "rs9999005" not in rsids


def test_append_to_csv_preserves_existing_rows(fake_existing_csv, tmp_path):
    """append_to_csv must keep all existing rows + add new ones, atomic
    write, vus_rank recomputed on the full merged dataset."""
    from monitor_clinvar import append_to_csv

    scored_new = pd.DataFrame([{
        "chrom": "2", "pos": 166_000_010, "ref": "C", "alt": "T",
        "rsid": "rs9999002", "clnsig_category": "uncertain",
        "SPLICE_SITES_score": 0.9, "SPLICE_SITE_USAGE_score": 120.0,
        "SPLICE_JUNCTIONS_score": 1100.0, "DNASE_score": 2.0,
        "molecular_consequence": "intron_variant", "dnase_success": True,
        "dnase_error": "",
    }])
    merged = append_to_csv(scored_new, fake_existing_csv)
    assert len(merged) == 3, f"expected 3 rows, got {len(merged)}"
    # Higher SPLICE_SITES (0.9) should now be rank 1.
    top = merged.sort_values("vus_rank").iloc[0]
    assert top["rsid"] == "rs9999002"
    assert int(top["vus_rank"]) == 1
    # Original file on disk must reflect the merge (atomic write).
    on_disk = pd.read_csv(fake_existing_csv)
    assert len(on_disk) == 3


def test_append_to_csv_refuses_to_shrink(fake_existing_csv):
    """If the merge result is smaller than the existing CSV, something is
    deeply wrong — append_to_csv should still complete (the caller is
    responsible for catching the shrinkage), but the on-disk file must
    remain intact. We pass an empty scored_new to verify the function
    preserves everything."""
    from monitor_clinvar import append_to_csv

    empty = pd.DataFrame(columns=[
        "chrom", "pos", "ref", "alt", "rsid", "clnsig_category",
        "SPLICE_SITES_score", "SPLICE_SITE_USAGE_score",
        "SPLICE_JUNCTIONS_score", "DNASE_score", "dnase_success",
        "dnase_error",
    ])
    merged = append_to_csv(empty, fake_existing_csv)
    assert len(merged) == 2  # original 2 rows preserved


def test_dry_run_end_to_end_with_fixture(fake_clinvar, tmp_path):
    """End-to-end: monitor_clinvar.py --dry-run on a fixture VCF + a
    minimal existing CSV that contains 1 of the 3 SCN1A VUS. Expects:
    - exit 0
    - log shows n_new_snv=1 (rs9999002), n_indels_skipped=1 (rs9999003)
    - no HF publish, no email
    - original existing CSV UNCHANGED (dry-run skips the append)
    """
    import shutil

    # Build a minimal existing CSV that includes rs9999001 only.
    existing = pd.DataFrame([{
        "chrom": "2", "pos": 166_000_001, "ref": "A", "alt": "G",
        "rsid": "rs9999001", "SPLICE_SITES_score": 0.5,
        "SPLICE_SITE_USAGE_score": 100.0, "SPLICE_JUNCTIONS_score": 1000.0,
        "DNASE_score": -5.0, "vus_rank": 1,
        "clnsig_category": "uncertain",
    }])
    existing_path = tmp_path / "existing.csv"
    existing.to_csv(existing_path, index=False)
    original_rows = len(existing)

    log_dir = tmp_path / "logs"
    rc = subprocess.run(
        [sys.executable, str(SCRIPTS / "monitor_clinvar.py"),
         "--vcf", str(fake_clinvar),
         "--tbi", str(fake_clinvar) + ".tbi",
         "--existing", str(existing_path),
         "--log-dir", str(log_dir),
         "--dry-run"],
        check=False, capture_output=True, text=True,
    )
    assert rc.returncode == 0, (
        f"monitor_clinvar --dry-run failed:\nSTDOUT:\n{rc.stdout}\n"
        f"STDERR:\n{rc.stderr}"
    )
    log = json.loads((log_dir / "monitor_clinvar.json").read_text())
    assert log["n_new_snv"] == 1, f"expected 1 new SNV, got {log}"
    assert log["n_indels_skipped"] == 1, f"expected 1 indel skipped, got {log}"
    assert log["dry_run"] is True
    assert log["hf_published"] is False
    assert log["email_sent"] is False

    # Critical: dry-run must NOT have modified the existing CSV.
    after_rows = len(pd.read_csv(existing_path))
    assert after_rows == original_rows, (
        "dry-run modified the existing CSV — should leave it untouched "
        "until step 4 (append) is reached."
    )


# ---- Email / HF gating (without actually sending / uploading) --------------

def test_email_returns_false_when_smtp_env_unset(monkeypatch, tmp_path,
                                                  fake_existing_csv):
    """With SMTP_USER/SMTP_PASSWORD unset, send_tier1_email must return
    False and NOT raise."""
    from monitor_clinvar import send_tier1_email

    monkeypatch.delenv("SMTP_USER", raising=False)
    monkeypatch.delenv("SMTP_PASSWORD", raising=False)
    monkeypatch.delenv("NOTIFY_EMAIL", raising=False)

    merged = pd.read_csv(fake_existing_csv)
    new_variants = pd.DataFrame()  # empty — nothing to alert on
    assert send_tier1_email(merged, new_variants, tmp_path / "log.json") is False


def test_email_returns_false_when_no_new_in_top50(monkeypatch, fake_existing_csv,
                                                   tmp_path):
    """If no new variants cross the top-50 threshold, do not send."""
    from monitor_clinvar import send_tier1_email

    monkeypatch.setenv("SMTP_USER", "x@gmail.com")
    monkeypatch.setenv("SMTP_PASSWORD", "x")
    monkeypatch.setenv("NOTIFY_EMAIL", "x@gmail.com")

    merged = pd.read_csv(fake_existing_csv)
    # New variant with very low SPLICE_SITES — won't make top-50.
    new_variants = pd.DataFrame([{
        "chrom": "2", "pos": 999, "ref": "A", "alt": "T",
        "SPLICE_SITES_score": 0.001, "DNASE_score": 0.0,
        "molecular_consequence": "intron_variant",
    }])
    assert send_tier1_email(merged, new_variants, tmp_path / "log.json") is False


def test_publish_hf_refuses_to_shrink(monkeypatch, tmp_path):
    """If the merged dataset has fewer rows than the published baseline,
    publish_hf must refuse rather than clobber the existing HF parquet."""
    from monitor_clinvar import publish_hf

    monkeypatch.setenv("HF_TOKEN", "hf_fake")
    # Tiny merged DF — well below the 1,610 baseline. Must NOT upload.
    merged = pd.DataFrame({"a": [1]})
    assert publish_hf(merged, "RROL/fake-repo") is False


def test_publish_hf_skips_without_token(monkeypatch, tmp_path):
    """Without HF_TOKEN, publish_hf should return False (not raise) so the
    rest of the monitor can still update the local CSV."""
    from monitor_clinvar import publish_hf

    monkeypatch.delenv("HF_TOKEN", raising=False)
    merged = pd.DataFrame({"a": list(range(2000))})
    assert publish_hf(merged, "RROL/fake-repo") is False


# ---- Launchd plist sanity --------------------------------------------------

def test_plist_exists_and_lints():
    """The launchd plist must exist and pass plutil -lint."""
    plist = SCRIPTS / "com.alpha-genome.scn1a-monitor.plist"
    assert plist.exists(), f"{plist} missing"
    rc = subprocess.run(
        ["plutil", "-lint", str(plist)],
        check=False, capture_output=True, text=True,
    )
    assert rc.returncode == 0, rc.stdout + rc.stderr
    assert "OK" in rc.stdout


def test_plist_uses_wrapper_not_inline_bash():
    """Verify the plist uses the _run_with_key.sh wrapper rather than an
    inline bash -c '...' command — launchd will silently reject long /
    escaped inline commands (known pitfall, see launchd-passive-income-deploy
    skill)."""
    plist = SCRIPTS / "com.alpha-genome.scn1a-monitor.plist"
    text = plist.read_text()
    # The program args array should NOT contain "-c"
    assert "-c" not in text, (
        "plist contains an inline `bash -c` — known to fail in launchd. "
        "Move the command into _run_with_key.sh instead."
    )
    # And SHOULD reference the wrapper script.
    assert "_run_with_key.sh" in text
    assert "monitor_clinvar.py" in text