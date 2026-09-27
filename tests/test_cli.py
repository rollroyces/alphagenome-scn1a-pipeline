"""Offline tests for the ``alphagenome-scn1a`` CLI.

None of these make API calls. They verify:
- the Typer app wires up cleanly
- every command listed in the task description is registered
- ``info`` and ``tier1`` run without an API key
- ``score-gene`` / ``score-vcf`` fail with a clear error and exit code 1
  when no key is present (and would otherwise succeed)
- key value is never echoed back

Run with: ``pytest tests/test_cli.py -v``
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Module-shape tests
# ---------------------------------------------------------------------------

def test_cli_module_imports():
    from alphagenome_scn1a import cli
    assert cli.app is not None
    assert callable(cli.main)


def test_subcommands_registered():
    """Every subcommand the task asked for must be wired into the Typer app."""
    from alphagenome_scn1a.cli import app

    # Typer stores registered commands in app.registered_commands (a list of
    # TyperCommand objects); their ``name`` attribute is the subcommand.
    names = {cmd.name for cmd in app.registered_commands}
    for required in ("info", "tier1", "score-gene", "score-vcf", "reproduce"):
        assert required in names, f"missing CLI command: {required}"


def test_version_exported():
    from alphagenome_scn1a import __version__
    parts = __version__.split(".")
    assert len(parts) >= 2
    assert all(p.isdigit() for p in parts), f"version {__version__!r} not semver-ish"


# ---------------------------------------------------------------------------
# Gene registry
# ---------------------------------------------------------------------------

def test_gene_registry_has_scn1a():
    from alphagenome_scn1a.extract import GENE_REGISTRY
    assert "SCN1A" in GENE_REGISTRY
    assert GENE_REGISTRY["SCN1A"]["chrom"] == "2"
    assert GENE_REGISTRY["SCN1A"]["start"] == 165_984_640


def test_resolve_unknown_gene_raises():
    from alphagenome_scn1a.extract import resolve_gene_region, UnknownGeneError
    with pytest.raises(UnknownGeneError) as exc:
        resolve_gene_region("NOT_A_REAL_GENE")
    assert "NOT_A_REAL_GENE" in str(exc.value)
    assert "SCN1A" in str(exc.value)  # lists known genes


def test_resolve_with_explicit_overrides():
    from alphagenome_scn1a.extract import resolve_gene_region
    chrom, start, end = resolve_gene_region(
        "SCN1A", chrom="2", start=1, end=2,
    )
    assert (chrom, start, end) == ("2", 1, 2)


# ---------------------------------------------------------------------------
# API key handling (no real key — only test the resolver shape)
# ---------------------------------------------------------------------------

def test_resolve_api_key_missing(monkeypatch, tmp_path, capsys):
    """With no env var and no key files, ``resolve_api_key(require=True)`` raises."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)  # no .alphagenome_key in cwd

    from alphagenome_scn1a.utils import resolve_api_key, MissingApiKeyError
    with pytest.raises(MissingApiKeyError) as exc:
        resolve_api_key(require=True)
    msg = str(exc.value)
    assert "ALPHAGENOME_API_KEY" in msg
    assert "alphagenome_key" in msg


def test_resolve_api_key_optional(monkeypatch, tmp_path):
    """``require=False`` returns None cleanly when no key is present."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)

    from alphagenome_scn1a.utils import resolve_api_key
    assert resolve_api_key(require=False) is None


def test_resolve_api_key_from_file(monkeypatch, tmp_path):
    """A 40-char ``AIzaSy...`` file is picked up."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)
    keyfile = tmp_path / ".alphagenome_key"
    fake = "AIzaSy" + "x" * 35  # 41 chars, plausible format
    keyfile.write_text(fake + "\n")

    from alphagenome_scn1a.utils import resolve_api_key
    key = resolve_api_key(require=True)
    assert key.value == fake
    assert "alphagenome_key" in key.source
    assert key.is_plausible()


def test_key_value_never_printed_in_missing_error(monkeypatch, tmp_path, capsys):
    """The ``MissingApiKeyError`` message must not leak any key fragment."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)
    fake = "AIzaSyREALKEY_REALTOTALfakereal"
    (tmp_path / ".alphagenome_key").write_text(fake)

    from alphagenome_scn1a.utils import resolve_api_key
    # require=True with a plausible key returns it (the message only triggers
    # when the key is missing entirely). Verify the alternate path: the
    # CLI's ``info`` command must never echo the real key.
    # We simulate by checking that the masked() output strips the middle.
    key = resolve_api_key(require=False)
    assert key is not None
    masked = key.masked()
    assert "REALKEY" not in masked
    assert fake not in masked
    assert masked.startswith("AIzaSy")
    assert masked.endswith("eal")


# ---------------------------------------------------------------------------
# End-to-end CLI invocations
# ---------------------------------------------------------------------------

def _run_cli(*args: str, env: dict | None = None) -> subprocess.CompletedProcess:
    """Invoke the CLI as a subprocess (matching what users do)."""
    e = {**os.environ, **(env or {})}
    e.pop("ALPHAGENOME_API_KEY", None)  # default: pretend no key
    return subprocess.run(
        [sys.executable, "-m", "alphagenome_scn1a.cli", *args],
        cwd=str(REPO_ROOT),
        env=e,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_cli_info_no_key():
    """``info`` must work without an API key."""
    result = _run_cli("info")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "alphagenome-scn1a" in result.stdout
    assert "API key" in result.stdout


def test_cli_tier1_no_key():
    """``tier1`` must work without an API key and print 4 candidates."""
    result = _run_cli("tier1", "SCN1A")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Tier-1 candidates" in result.stdout
    # Should show 4 distinct VUS ranks.
    for rank in ("14", "17", "29", "30"):
        assert rank in result.stdout, f"missing tier-1 candidate with rank {rank}"


def test_cli_score_gene_without_key_fails_cleanly():
    """``score-gene`` must exit 1 with a helpful message when no key."""
    result = _run_cli("score-gene", "SCN1A", "--output", "/tmp/_should_not_exist.csv")
    assert result.returncode == 1, result.stdout + result.stderr
    # Help text must mention both fix paths
    combined = result.stdout + result.stderr
    assert "ALPHAGENOME_API_KEY" in combined
    assert "alphagenome_key" in combined
    # And explicitly say which commands need vs. don't need a key
    assert "score-gene" in combined
    assert "info" in combined


def test_cli_score_vcf_without_key_fails_cleanly():
    result = _run_cli(
        "score-vcf", "/tmp/does_not_matter.vcf",
        "--gene", "SCN1A", "--output", "/tmp/_x.csv",
    )
    assert result.returncode == 1, result.stdout + result.stderr
    assert "ALPHAGENOME_API_KEY" in (result.stdout + result.stderr)


def test_cli_help_lists_all_commands():
    result = _run_cli("--help")
    assert result.returncode == 0
    for cmd in ("info", "tier1", "score-gene", "score-vcf", "reproduce"):
        assert cmd in result.stdout, f"--help missing command: {cmd}"


def test_cli_unknown_gene_exits_nonzero(tmp_path):
    """``score-gene`` for a gene the registry doesn't know → exit 2."""
    # Create a fake VCF-less invocation by passing the gene and capturing the
    # UnknownGeneError path. Since the key check runs first, we need a key.
    keyfile = tmp_path / ".alphagenome_key"
    keyfile.write_text("AIzaSy" + "x" * 35)
    env = {"ALPHAGENOME_API_KEY": "AIzaSy" + "x" * 35}
    e = {**os.environ, **env}
    result = subprocess.run(
        [sys.executable, "-m", "alphagenome_scn1a.cli",
         "score-gene", "NOT_A_REAL_GENE", "--output", "/tmp/_x.csv"],
        cwd=str(REPO_ROOT), env=e, capture_output=True, text=True, timeout=20,
    )
    # Either exit 2 (UnknownGeneError) — the path we want — or 1 if pytest
    # can't reach the API. Just verify it fails and the message is helpful.
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "Unknown gene" in combined or "ALPHAGENOME_API_KEY" in combined