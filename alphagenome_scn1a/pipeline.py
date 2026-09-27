"""End-to-end pipeline: extract ClinVar variants → score → AUPRC → write CSV.

This is the high-level orchestrator that ``score-gene`` calls. It's a thin
sequence of calls into :mod:`alphagenome_scn1a.extract` and
:mod:`alphagenome_scn1a.score`; we keep it separate so the CLI handler stays
declarative and easy to skim.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import extract, score


def run_score_gene(
    gene: str,
    output_csv: str | Path,
    api_key: str,
    *,
    n_pos: int = 200,
    n_neg: int = 350,
    extract_tsv: str | Path | None = None,
    skip_extract: bool = False,
    on_progress=None,
) -> dict:
    """Run the full extract → score → metrics pipeline for one gene.

    Args:
        gene: gene symbol (e.g. ``"SCN1A"``). Coordinates looked up in
              :data:`extract.GENE_REGISTRY`.
        output_csv: where to write the final scored CSV (also contains the
                    AUPRC summary in a separate file alongside it).
        api_key: AlphaGenome API key.
        n_pos / n_neg: positive / negative cap during ClinVar sampling.
        extract_tsv: where to write (or read) the intermediate TSV. Defaults
                     to ``outputs/clinvar_{gene}_benchmark.tsv``.
        skip_extract: reuse the existing ``extract_tsv`` instead of rebuilding.
        on_progress: optional progress callback, see ``score.score_dataframe``.

    Returns:
        Summary dict with ``n_total``, ``n_success``, ``n_fail``, AUPRC values,
        and the paths to the written files. The CLI renders this to the user.
    """
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if extract_tsv is None:
        extract_tsv = output_csv.parent / f"clinvar_{gene.lower()}_benchmark.tsv"
    extract_tsv = Path(extract_tsv)

    # ---- Extract ---------------------------------------------------------
    n_variants = 0
    if not skip_extract or not extract_tsv.exists():
        n_variants = extract.extract_benchmark(
            gene, extract_tsv, n_pos=n_pos, n_neg=n_neg,
        )
    else:
        # Reuse the pre-extracted TSV.
        existing = pd.read_csv(extract_tsv, sep="\t")
        n_variants = len(existing)

    # ---- Score -----------------------------------------------------------
    summary = score.score_benchmark_file(
        extract_tsv, output_csv, api_key, gene=gene, on_progress=on_progress,
    )

    # ---- Metrics (AUPRC + AUROC + top-5% precision) -----------------------
    metrics = compute_summary_metrics(output_csv)
    summary.update({
        "extract_tsv": str(extract_tsv),
        "n_variants_extracted": n_variants,
        "metrics": metrics,
    })
    return summary


def compute_summary_metrics(scored_csv: str | Path) -> dict:
    """Compute AUROC, AUPRC, and top-5% precision for SPLICE_SITES scorer.

    Returns an empty dict if the file has no successful scores or only one
    class. The CLI renders "—" instead of crashing in that case.
    """
    df = pd.read_csv(scored_csv)
    if "score_success" not in df.columns or "label" not in df.columns:
        return {}
    ok = df[df["score_success"] == True].copy()  # noqa: E712 — pandas-idiomatic
    if ok.empty or ok["label"].nunique() < 2:
        return {"reason": "insufficient data (need successes from both labels)"}

    y_true = (ok["label"] == "positive").astype(int).values

    out: dict = {
        "n_total": int(len(ok)),
        "n_positive": int(y_true.sum()),
        "n_negative": int((1 - y_true).sum()),
    }
    if "SPLICE_SITES_score" in ok.columns:
        y_score = ok["SPLICE_SITES_score"].values
        try:
            from sklearn.metrics import average_precision_score, roc_auc_score
            out["SPLICE_SITES_auroc"] = float(roc_auc_score(y_true, y_score))
            out["SPLICE_SITES_auprc"] = float(average_precision_score(y_true, y_score))
        except Exception as e:  # pragma: no cover — sklearn missing in slim envs
            out["metrics_error"] = f"sklearn not available: {e}"

        # Top-5% precision
        import numpy as np
        n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
        top_idx = np.argsort(y_score)[-n_top:]
        out["top_5_pct_precision"] = float(y_true[top_idx].mean())

    return out


def render_metrics(metrics: dict) -> str:
    """Format the metrics dict as a multi-line, human-readable string."""
    if not metrics:
        return "  (no metrics available)"
    if "reason" in metrics:
        return f"  ({metrics['reason']})"
    lines = []
    if "n_total" in metrics:
        lines.append(
            f"  variants scored: {metrics['n_total']} "
            f"({metrics.get('n_positive', '?')} pos, {metrics.get('n_negative', '?')} neg)"
        )
    for key in ("SPLICE_SITES_auroc", "SPLICE_SITES_auprc", "top_5_pct_precision"):
        if key in metrics:
            label = key.replace("_", " ").replace("SPLICE SITES", "SPLICE_SITES")
            lines.append(f"  {label}: {metrics[key]:.4f}")
    return "\n".join(lines)