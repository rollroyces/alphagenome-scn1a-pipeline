"""Console entry point: ``alphagenome-scn1a <subcommand>``.

A friendly wrapper over the per-experiment scripts in ``scripts/`` and the
modules in :mod:`alphagenome_scn1a`. Built on Typer for self-documenting
help and clean error messages.

Subcommands
-----------
score-gene GENE --output PATH
    Extract ClinVar variants for ``GENE``, score them with AlphaGenome,
    write a scored CSV. AUPRC / AUROC are printed at the end.

score-vcf VCF --gene GENE --output PATH
    Score variants from a user-supplied VCF instead of pulling from ClinVar.

reproduce
    Run ``make reproduce`` — the offline + small-API-call re-verification
    harness.

info
    Show installed version, API key status. No API call.

tier1 [GENE]
    Print the cached Tier-1 candidates from ``outputs/tier1_actionable_features.csv``.
    No API call.

smoke / check (legacy)
    The original 0.1 shim — kept so existing Makefile / Docker recipes keep
    working until they're updated.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import typer
from typing_extensions import Annotated

from . import __version__
from .pipeline import run_score_gene, render_metrics
from .utils import (
    MissingApiKeyError,
    api_key_status_lines,
    banner,
    humanize_path,
    resolve_api_key,
)
from . import extract as _extract_mod


REPO_ROOT = Path(__file__).resolve().parent.parent

app = typer.Typer(
    name="alphagenome-scn1a",
    help=(
        "AlphaGenome SCN1A / rare-disease variant scoring CLI.\n\n"
        "Quick examples:\n"
        "  alphagenome-scn1a info\n"
        "  alphagenome-scn1a tier1 SCN1A\n"
        "  alphagenome-scn1a score-gene SCN1A --output outputs/scored_scn1a.csv"
    ),
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode=None,
)


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@app.command()
def info() -> None:
    """Show installed version, SDK presence, and API key status.

    Does NOT make any API calls — safe to run before you have a key.
    """
    typer.echo(banner("alphagenome-scn1a"))
    typer.echo(f"  version:        {__version__}")
    typer.echo(f"  python:         {sys.version.split()[0]}")
    typer.echo(f"  repo root:      {humanize_path(REPO_ROOT)}")
    typer.echo()

    # SDK presence
    try:
        import alphagenome  # noqa: F401
        typer.echo("  alphagenome:    installed")
    except ImportError as e:
        typer.echo(f"  alphagenome:    MISSING ({e})")

    # Typer presence
    typer.echo(f"  typer:          installed")
    typer.echo()

    typer.echo(banner("API key"))
    for line in api_key_status_lines():
        typer.echo(f"  {line}")
    typer.echo()

    typer.echo(banner("Quickstart"))
    typer.echo("  alphagenome-scn1a tier1 SCN1A     # no API call")
    typer.echo("  alphagenome-scn1a score-gene SCN1A --output results.csv")
    typer.echo("  alphagenome-scn1a reproduce       # verify existing outputs")


# ---------------------------------------------------------------------------
# tier1
# ---------------------------------------------------------------------------

@app.command()
def tier1(
    gene: Annotated[str, typer.Argument(help="Gene symbol, e.g. SCN1A.")] = "SCN1A",
) -> None:
    """Print the cached Tier-1 candidates for ``GENE`` (no API call).

    Reads from ``outputs/tier1_actionable_features.csv`` which the paper
    calls the "4 Tier-1 candidates ready for outreach" — splice-altering
    variants that match multiple pathogenicity criteria.
    """
    csv_path = REPO_ROOT / "outputs" / "tier1_actionable_features.csv"
    if not csv_path.exists():
        typer.echo(f"ERROR: {csv_path} not found.", err=True)
        typer.echo(
            "Run the full pipeline (make reproduce) to generate it.",
            err=True,
        )
        raise typer.Exit(code=1)

    import pandas as pd
    df = pd.read_csv(csv_path)
    if df.empty:
        typer.echo(f"No Tier-1 candidates in {csv_path.name}.", err=True)
        raise typer.Exit(code=1)

    typer.echo(banner(f"Tier-1 candidates for {gene}"))
    typer.echo(f"  source:  {humanize_path(csv_path)}")
    typer.echo(f"  count:   {len(df)}")
    typer.echo()

    # Show a tight summary table — don't dump every column.
    cols = [
        "vus_rank", "rsid", "chrom", "pos", "ref", "alt",
        "molecular_consequence", "clndn", "SPLICE_SITES_score",
    ]
    cols = [c for c in cols if c in df.columns]
    rows = df[cols].to_string(index=False)
    typer.echo(rows)
    typer.echo()
    typer.echo(
        f"Full table (with DNASE / gnomAD / RNA-seq): {humanize_path(csv_path)}"
    )


# ---------------------------------------------------------------------------
# score-gene
# ---------------------------------------------------------------------------

@app.command()
def score_gene(
    gene: Annotated[str, typer.Argument(help="Gene symbol (e.g. SCN1A, KCNQ2, COL4A5, FBN1).")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Where to write the scored CSV.")] = Path("outputs/scored.csv"),
    n_pos: Annotated[int, typer.Option(help="Max positive (pathogenic) controls.")] = 200,
    n_neg: Annotated[int, typer.Option(help="Max negative (benign) controls.")] = 350,
    skip_extract: Annotated[bool, typer.Option(
        "--skip-extract",
        help="Reuse a pre-extracted TSV at outputs/clinvar_<gene>_benchmark.tsv.",
    )] = False,
) -> None:
    """Extract ClinVar variants for ``GENE``, score with AlphaGenome, write CSV.

    Requires an API key (set ``$ALPHAGENOME_API_KEY`` or write
    ``~/.alphagenome_key``).
    """
    try:
        key = resolve_api_key(require=True)
    except MissingApiKeyError as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1)

    if not output.is_absolute():
        output = REPO_ROOT / output

    typer.echo(banner(f"score-gene {gene.upper()}"))
    typer.echo(f"  output:       {humanize_path(output)}")
    typer.echo(f"  api key:      {key.masked()}  (source: {key.source})")
    typer.echo(f"  extract caps: n_pos={n_pos}, n_neg={n_neg}")
    typer.echo()

    def _on_progress(i, n, n_success, n_fail, elapsed, eta) -> None:
        typer.echo(
            f"  [{i+1}/{n}] {elapsed:.0f}s elapsed, "
            f"success={n_success}, fail={n_fail}, ETA {eta:.0f}s"
        )

    try:
        summary = run_score_gene(
            gene.upper(),
            output,
            key.value,
            n_pos=n_pos,
            n_neg=n_neg,
            skip_extract=skip_extract,
            on_progress=_on_progress,
        )
    except _extract_mod.UnknownGeneError as e:
        typer.echo(f"ERROR: {e}", err=True)
        raise typer.Exit(code=2)

    typer.echo()
    typer.echo(banner("Result"))
    typer.echo(f"  variants extracted: {summary['n_variants_extracted']}")
    typer.echo(f"  variants scored:    {summary['n_success']}/{summary['n_total']}")
    if summary["n_fail"]:
        typer.echo(f"  variants failed:    {summary['n_fail']}")
    typer.echo(f"  raw CSV:            {humanize_path(summary['output'])}")
    typer.echo()
    typer.echo(banner("Metrics"))
    typer.echo(render_metrics(summary["metrics"]))


# ---------------------------------------------------------------------------
# score-vcf
# ---------------------------------------------------------------------------

@app.command()
def score_vcf(
    vcf: Annotated[Path, typer.Argument(help="Path to a VCF file (bgzipped or plain).")],
    gene: Annotated[str, typer.Option("--gene", help="Gene symbol, written into the output rows.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Where to write the scored CSV.")] = Path("outputs/scored_vcf.csv"),
) -> None:
    """Score variants from a custom VCF with AlphaGenome.

    The VCF must have ``CHROM POS REF ALT`` columns (standard VCF). Multi-allelic
    records are split into one row per ALT allele.
    """
    try:
        key = resolve_api_key(require=True)
    except MissingApiKeyError as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1)

    if not vcf.exists():
        typer.echo(f"ERROR: VCF not found: {vcf}", err=True)
        raise typer.Exit(code=2)
    if not output.is_absolute():
        output = REPO_ROOT / output

    typer.echo(banner(f"score-vcf {vcf.name}"))
    typer.echo(f"  gene:        {gene.upper()}")
    typer.echo(f"  output:      {humanize_path(output)}")
    typer.echo(f"  api key:     {key.masked()}  (source: {key.source})")
    typer.echo()

    variants = _read_vcf(vcf)
    if not variants:
        typer.echo("ERROR: No SNV variants found in VCF.", err=True)
        raise typer.Exit(code=3)

    typer.echo(f"  loaded {len(variants)} SNV variants")

    from . import score as _score_mod

    def _on_progress(i, n, n_success, n_fail, elapsed, eta) -> None:
        typer.echo(
            f"  [{i+1}/{n}] {elapsed:.0f}s elapsed, "
            f"success={n_success}, fail={n_fail}, ETA {eta:.0f}s"
        )

    df = _score_mod.score_variants(
        variants, key.value, gene=gene.upper(), on_progress=_on_progress,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    n_success = int((df["score_success"] == True).sum()) if not df.empty else 0  # noqa: E712
    typer.echo()
    typer.echo(banner("Result"))
    typer.echo(f"  scored:      {n_success}/{len(df)}")
    typer.echo(f"  raw CSV:     {humanize_path(output)}")


def _read_vcf(path: Path) -> list[dict]:
    """Parse a (possibly bgzipped) VCF into one dict per SNV allele.

    Skips header lines (``#``) and indels. Doesn't require pysam so a plain
    VCF works the same as a ``.vcf.gz``.
    """
    opener = _open_maybe_gzip(path)
    out: list[dict] = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _id, ref, alts = fields[:5]
            for alt in alts.split(","):
                if len(ref) == 1 and len(alt) == 1 and alt != "*":
                    out.append({
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                    })
    return out


def _open_maybe_gzip(path: Path):
    """Return a callable that opens ``path`` as text, gzip-decompressing if needed."""
    if str(path).endswith(".gz"):
        import gzip
        return lambda p: gzip.open(p, "rt")
    return lambda p: open(p, "r")


# ---------------------------------------------------------------------------
# reproduce
# ---------------------------------------------------------------------------

@app.command()
def reproduce() -> None:
    """Run ``make reproduce`` — re-verify all outputs without burning API quota.

    This delegates to the Makefile target defined by the reproducibility
    harness subagent. Fails fast if any expected output is missing or if
    any small re-run (smoke, ISM) fails.
    """
    typer.echo(banner("reproduce"))
    typer.echo("Delegating to `make reproduce`...")
    typer.echo()
    rc = subprocess.call(["make", "reproduce"], cwd=str(REPO_ROOT))
    if rc != 0:
        typer.echo()
        typer.echo("reproduce FAILED — see make output above.", err=True)
        raise typer.Exit(code=rc)
    typer.echo()
    typer.echo("reproduce PASSED.")


# ---------------------------------------------------------------------------
# legacy: smoke / check (kept so old Makefile / Docker recipes work)
# ---------------------------------------------------------------------------

@app.command(hidden=True)
def smoke() -> None:
    """Run scripts/smoke_test.py (legacy alias)."""
    cmd = [sys.executable, str(REPO_ROOT / "scripts" / "smoke_test.py")]
    typer.echo(f"[alphagenome-scn1a] {' '.join(cmd)}")
    rc = subprocess.call(cmd, cwd=str(REPO_ROOT))
    raise typer.Exit(code=rc)


@app.command(hidden=True)
def check() -> None:
    """Run ``make check`` (legacy alias)."""
    rc = subprocess.call(["make", "check"], cwd=str(REPO_ROOT))
    raise typer.Exit(code=rc)


# ---------------------------------------------------------------------------
# Typer entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Console-script entry point declared in pyproject.toml."""
    try:
        app()
    except KeyboardInterrupt:
        typer.echo("\nInterrupted.", err=True)
        raise typer.Exit(code=130) from None


if __name__ == "__main__":
    main()