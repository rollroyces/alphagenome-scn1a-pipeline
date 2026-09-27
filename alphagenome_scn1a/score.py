"""Wrap ``scripts/score_gene_benchmark.py`` so the CLI can call it.

We re-use the existing ``score_one`` / bootstrapping logic rather than
duplicating it. The wrapper exposes only the operations the CLI needs:
load a TSV, score it with AlphaGenome, write the raw CSV.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd


SCORER_NAMES = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]


def score_dataframe(
    df: pd.DataFrame,
    api_key: str,
    *,
    gene: str | None = None,
    on_progress: Callable[[int, int, int, int, float, float], None] | None = None,
) -> pd.DataFrame:
    """Score every row of ``df`` with AlphaGenome.

    ``df`` must have columns: chrom, pos, ref, alt. Optional columns that
    are carried through: rsid, label, clnsig_category, molecular_consequence.

    Args:
        df: input benchmark (one row per variant).
        api_key: AlphaGenome API key.
        gene: overrides the ``gene`` column in the output. Defaults to the
              first row's ``gene`` value, or "unknown".
        on_progress: optional ``callable(i, n, n_success, n_fail, elapsed, eta)``
              invoked every 25 rows so the CLI can print a progress line.

    Returns:
        A new DataFrame with the original columns plus SPLICE_*_score,
        score_success, and (on failure) score_error.
    """
    # Import inside the function so `info` / `tier1` don't require alphagenome.
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    dna_model = dna_client.create(api_key)
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[n] for n in SCORER_NAMES]

    gene_label = gene or (df["gene"].iloc[0] if "gene" in df.columns and len(df) else "unknown")
    n = len(df)
    rows: list[dict] = []
    t0 = time.time()
    n_success = 0
    n_fail = 0

    for i, row in df.iterrows():
        if i % 25 == 0 or i == n - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0.0
            eta = (n - i - 1) / rate if rate > 0 else 0.0
            if on_progress:
                on_progress(i, n, n_success, n_fail, elapsed, eta)

        result = _score_one(
            dna_model, scorers,
            row["chrom"], row["pos"], row["ref"], row["alt"],
        )
        result.update({
            "gene": gene_label,
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "rsid": row.get("rsid", "") if hasattr(row, "get") else row.get("rsid", ""),
            "label": row.get("label", "") if hasattr(row, "get") else row.get("label", ""),
            "clnsig_category": row.get("clnsig_category", "") if hasattr(row, "get") else row.get("clnsig_category", ""),
            "molecular_consequence": row.get("molecular_consequence", "") if hasattr(row, "get") else row.get("molecular_consequence", ""),
        })
        rows.append(result)
        if result["score_success"]:
            n_success += 1
        else:
            n_fail += 1

    return pd.DataFrame(rows)


def _score_one(dna_model, scorers, chrom, pos, ref, alt) -> dict:
    """Score a single variant; on failure return a row with score_success=False."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    try:
        chrom_str = str(chrom)
        if not chrom_str.startswith("chr"):
            chrom_str = "chr" + chrom_str
        variant = genome.Variant(
            chromosome=chrom_str,
            position=int(pos),
            reference_bases=str(ref),
            alternate_bases=str(alt),
        )
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
        scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
        )
        out: dict = {"score_success": True}
        for ann, name in zip(scores, SCORER_NAMES):
            out[f"{name}_score"] = float(ann.X.sum())
        return out
    except Exception as e:
        return {
            "score_success": False,
            "score_error": f"{type(e).__name__}: {str(e)[:120]}",
        }


def score_benchmark_file(
    input_tsv: str | Path,
    output_csv: str | Path,
    api_key: str,
    *,
    gene: str | None = None,
    on_progress: callable | None = None,
) -> dict:
    """Score a TSV file end-to-end. Returns a small summary dict."""
    input_tsv = Path(input_tsv)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_tsv, sep="\t")
    scored = score_dataframe(df, api_key, gene=gene, on_progress=on_progress)
    scored.to_csv(output_csv, index=False)

    n = len(scored)
    n_success = int((scored["score_success"] == True).sum())
    return {
        "input": str(input_tsv),
        "output": str(output_csv),
        "n_total": n,
        "n_success": n_success,
        "n_fail": n - n_success,
    }


def score_variants(
    variants: Iterable[dict],
    api_key: str,
    *,
    gene: str = "unknown",
    on_progress: callable | None = None,
) -> pd.DataFrame:
    """Score an arbitrary list of variants. Each dict needs chrom/pos/ref/alt.

    Used by ``score-vcf`` after parsing a VCF.
    """
    df = pd.DataFrame(list(variants))
    if df.empty:
        return df
    return score_dataframe(df, api_key, gene=gene, on_progress=on_progress)