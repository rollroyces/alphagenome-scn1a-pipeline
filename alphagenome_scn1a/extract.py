"""Wrap ``scripts/extract_clinvar_for_gene.py`` so the CLI can call it.

The CLI only needs one knob — the gene name. We keep a small registry of
gene → (chrom, start, end) here; for unknown genes the user must pass
``--chrom``, ``--start``, ``--end`` explicitly.

This module deliberately does NOT duplicate the heavy lifting in the
scripts/ files. It re-imports ``extract_clinvar_for_gene.make_benchmark``
and ``fetch_variants`` so we share one source of truth with the per-experiment
scripts the researchers already trust.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make scripts/ importable so we can reuse the existing extraction logic.
REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPTS = REPO_ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


# GRCh38 coordinates for the genes covered in this repo. Source: GENCODE v46
# / NCBI RefSeq, all consistent with the values used in scripts/_*_extract.py.
#
# Each entry is the canonical window for AlphaGenome's 16 KB reference interval;
# flanking 50 KB on each side is plenty since the SDK resizes to 16 KB internally.
GENE_REGISTRY: dict[str, dict[str, int | str]] = {
    "SCN1A": {"chrom": "2",  "start": 165_984_640, "end": 166_182_806},
    "SCN2A": {"chrom": "2",  "start": 165_295_770, "end": 165_404_815},
    "KCNQ2": {"chrom": "20", "start":   63_400_679, "end":   63_472_909},
    "COL4A5": {"chrom": "X", "start":  108_439_837, "end":  108_697_545},
    "FBN1":   {"chrom": "15", "start":  48_408_312, "end":   48_645_721},
    "CFTR":   {"chrom": "7",  "start":  117_287_120, "end":  117_715_971},
    "DMD":    {"chrom": "X",  "start":   31_096_780, "end":   33_629_711},
    "MECP2":  {"chrom": "X",  "start":  154_021_573, "end":  154_097_717},
}


def resolve_gene_region(
    gene: str,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> tuple[str, int, int]:
    """Return (chrom, start, end) for ``gene``.

    Looks up ``gene`` in :data:`GENE_REGISTRY` first; if any of ``chrom``,
    ``start``, ``end`` are passed explicitly they override the registry values
    (useful for non-canonical transcripts or liftover coordinates).
    """
    gene_upper = gene.upper()
    if gene_upper not in GENE_REGISTRY and (chrom is None or start is None or end is None):
        raise UnknownGeneError(gene)
    base = GENE_REGISTRY.get(gene_upper, {})
    return (
        str(chrom if chrom is not None else base.get("chrom", "")),
        int(start if start is not None else base["start"]),
        int(end if end is not None else base["end"]),
    )


class UnknownGeneError(ValueError):
    def __init__(self, gene: str) -> None:
        self.gene = gene
        known = ", ".join(sorted(GENE_REGISTRY))
        super().__init__(
            f"Unknown gene: {gene!r}. Known genes: {known}. "
            f"Or pass --chrom / --start / --end explicitly."
        )


def extract_benchmark(
    gene: str,
    output_tsv: str | Path,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    n_pos: int = 200,
    n_neg: int = 350,
    seed: int = 42,
) -> int:
    """Extract a ClinVar benchmark for ``gene`` and write it to ``output_tsv``.

    Returns the number of variants written. Wraps the existing
    ``extract_clinvar_for_gene.make_benchmark`` so the per-experiment scripts
    and the CLI agree on the sampling logic.

    The output TSV has the same columns as the per-gene benchmark scripts:
    gene, chrom, pos, rsid, ref, alt, label, clnsig_category, molecular_consequence.
    """
    # Import lazily so missing pysam doesn't break `info` / `tier1`.
    import extract_clinvar_for_gene  # type: ignore[import-not-found]

    chrom_s, start_i, end_i = resolve_gene_region(gene, chrom, start, end)

    benchmark = extract_clinvar_for_gene.make_benchmark(
        chrom_s, start_i, end_i, gene.upper(),
        n_pos=n_pos, n_neg=n_neg, seed=seed,
    )
    if not benchmark:
        raise RuntimeError(f"No ClinVar variants found for {gene.upper()}")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    import csv
    with output_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "gene", "chrom", "pos", "rsid", "ref", "alt",
            "label", "clnsig_category", "molecular_consequence",
        ])
        for b in benchmark:
            v = b["variant"]
            w.writerow([
                gene.upper(), v["chrom"], v["pos"], v["rsid"],
                v["ref"], v["alts"][0], b["label"],
                v["clnsig_category"], v["molecular_consequence"],
            ])

    return len(benchmark)


def load_clinvar_table(gene: str) -> "Path | None":
    """Return the canonical pre-extracted ClinVar TSV for a gene, if it exists.

    The pipeline relies on already-extracted tables in ``outputs/`` so users
    don't pay the ClinVar query cost on every run. We look for the standard
    filenames used by the per-gene scripts.
    """
    candidates = [
        REPO_ROOT / "outputs" / f"clinvar_{gene.lower()}_benchmark.tsv",
        REPO_ROOT / "outputs" / f"clinvar_{gene.lower()}.tsv",
        REPO_ROOT / "outputs" / f"_{gene.lower()}_stratified.tsv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None