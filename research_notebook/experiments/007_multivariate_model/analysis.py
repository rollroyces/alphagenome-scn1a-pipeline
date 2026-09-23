"""
Experiment 007 — Multivariate model: DNase + ATAC + splicing → pathogenic vs benign.

Combines per-variant ISM metrics from Exp 001/002/003/004 and evaluates a
logistic regression with 5-fold stratified CV. Reports single-feature AUPRC,
multivariate AUPRC, naive baselines (random, SPLICE_SITES_score).

Outputs (all in this directory):
  - features.csv             : per-variant feature matrix with labels
  - cv_results.csv           : per-fold metrics for each model
  - comparison_table.md      : human-readable comparison
  - results.json             : machine-readable summary
  - README.md                : honest writeup
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=ConvergenceWarning)

REPO = Path("/Users/hermes/projects/alphagenome-work")
OUT = REPO / "research_notebook/experiments/007_multivariate_model"
OUT.mkdir(parents=True, exist_ok=True)

EXP001 = REPO / "research_notebook/experiments/001_ism_scn1a/metrics.csv"
EXP002 = REPO / "research_notebook/experiments/002_ism_atac/metrics.csv"
EXP003 = REPO / "research_notebook/experiments/003_ism_dnase_replication/metrics.csv"
EXP004 = REPO / "research_notebook/experiments/004_ism_dnase_n30/metrics.csv"

SPLICE_CSVS = {
    "scn1a": REPO / "outputs/benchmark_scn1a_live_api_raw.csv",
    "cftr":  REPO / "outputs/cross_disease_cftr_raw.csv",
    "dmd":   REPO / "outputs/cross_disease_dmd_raw.csv",
}

RANDOM_STATE = 42
N_SPLITS = 5


# ---------------------------------------------------------------------------
# 1. Load and merge ISM metrics from each experiment
# ---------------------------------------------------------------------------

def load_exp001() -> pd.DataFrame:
    """Splicing ISM (Exp 001): SCN1A only, n=10+10. Splice modality."""
    df = pd.read_csv(EXP001)
    df = df.rename(columns={
        "total":                    "splice_magnitude",
        "fraction_within_5bp":      "splice_frac_5bp",
        "fraction_within_15bp":     "splice_frac_15bp",
        "concentration_5bp":        "splice_concentration_5bp",
        "max_position_relative_to_variant": "splice_max_pos_rel",
        "max_value":                "splice_max_value",
    })
    df["gene"] = "scn1a"
    keep = ["gene", "chrom", "pos", "ref", "variant_label",
            "splice_magnitude", "splice_frac_5bp", "splice_frac_15bp",
            "splice_concentration_5bp", "splice_max_pos_rel", "splice_max_value"]
    return df[keep]


def load_exp002() -> pd.DataFrame:
    """DNase + ATAC ISM (Exp 002): SCN1A only, n=10+10, two modalities."""
    df = pd.read_csv(EXP002)
    df["gene"] = "scn1a"
    keep = ["gene", "chrom", "pos", "ref", "variant_label", "modality",
            "total", "fraction_within_5bp", "fraction_within_15bp",
            "concentration_5bp", "max_position_relative_to_variant", "max_value"]
    df = df[keep]

    dnase = (df[df["modality"] == "DNASE"]
             .drop(columns=["modality"])
             .rename(columns={
                 "total":                   "dnase_magnitude",
                 "fraction_within_5bp":     "dnase_frac_5bp",
                 "fraction_within_15bp":    "dnase_frac_15bp",
                 "concentration_5bp":       "dnase_concentration_5bp",
                 "max_position_relative_to_variant": "dnase_max_pos_rel",
                 "max_value":               "dnase_max_value",
             }))
    atac = (df[df["modality"] == "ATAC"]
            .drop(columns=["modality"])
            .rename(columns={
                "total":                   "atac_magnitude",
                "fraction_within_5bp":     "atac_frac_5bp",
                "fraction_within_15bp":    "atac_frac_15bp",
                "concentration_5bp":       "atac_concentration_5bp",
                "max_position_relative_to_variant": "atac_max_pos_rel",
                "max_value":               "atac_max_value",
            }))
    merged = dnase.merge(atac, on=["gene", "chrom", "pos", "ref", "variant_label"],
                         how="outer")
    return merged


def load_exp003() -> pd.DataFrame:
    """DNase ISM (Exp 003): DMD + CFTR, n=8+10 / 7+10."""
    df = pd.read_csv(EXP003)
    df = df.rename(columns={
        "total":                    "dnase_magnitude",
        "fraction_within_5bp":      "dnase_frac_5bp",
        "fraction_within_15bp":     "dnase_frac_15bp",
        "concentration_5bp":        "dnase_concentration_5bp",
        "max_position_relative_to_variant": "dnase_max_pos_rel",
        "max_value":                "dnase_max_value",
    })
    keep = ["gene", "chrom", "pos", "ref", "variant_label",
            "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
            "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value"]
    return df[keep]


def load_exp004() -> pd.DataFrame:
    """DNase ISM (Exp 004, primary): SCN1A + DMD + CFTR, n=30+30 (per gene)."""
    df = pd.read_csv(EXP004)
    df = df.rename(columns={
        "total":                    "dnase_magnitude",
        "fraction_within_5bp":      "dnase_frac_5bp",
        "fraction_within_15bp":     "dnase_frac_15bp",
        "concentration_5bp":        "dnase_concentration_5bp",
        "max_position_relative_to_variant": "dnase_max_pos_rel",
        "max_value":                "dnase_max_value",
    })
    keep = ["gene", "chrom", "pos", "ref", "variant_label",
            "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
            "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value"]
    return df[keep]


def _aggregate_alts(df: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Aggregate (mean) across alt alleles so we have one row per
    (gene, chrom, pos, ref, variant_label). Variants with the same SNV
    location but different alts appear as one row in the feature table.
    """
    return (df.groupby(["gene", "chrom", "pos", "ref", "variant_label"],
                       as_index=False, dropna=False)
              .agg({c: "mean" for c in value_cols}))


def build_feature_table() -> pd.DataFrame:
    # Aggregate Exp 001/002/004 to one row per (gene, chrom, pos, ref, label).
    # (Each of those files has multiple rows per SNV — one per alt.)
    e001 = _aggregate_alts(load_exp001(), [
        "splice_magnitude", "splice_frac_5bp", "splice_frac_15bp",
        "splice_concentration_5bp", "splice_max_pos_rel", "splice_max_value",
    ])
    e002 = _aggregate_alts(load_exp002(), [
        "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
        "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value",
        "atac_magnitude", "atac_frac_5bp", "atac_frac_15bp",
        "atac_concentration_5bp", "atac_max_pos_rel", "atac_max_value",
    ])
    e003 = _aggregate_alts(load_exp003(), [
        "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
        "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value",
    ])
    e004 = _aggregate_alts(load_exp004(), [
        "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
        "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value",
    ])

    # DNase: prefer Exp 004 (largest, all 3 genes), backfill with Exp 003 (DMD/CFTR),
    # then Exp 002 (SCN1A). Exp 004 should already cover everything.
    dnase_keys = ["gene", "chrom", "pos", "ref", "variant_label",
                  "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
                  "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value"]
    dnase = (e004[dnase_keys]
             .merge(e003[dnase_keys], on=["gene", "chrom", "pos", "ref", "variant_label"],
                    how="outer", suffixes=("", "_e3"))
             .merge(e002[["gene", "chrom", "pos", "ref", "variant_label",
                          "dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
                          "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value"]],
                    on=["gene", "chrom", "pos", "ref", "variant_label"],
                    how="outer", suffixes=("", "_e2")))
    # Coalesce e3 -> e4 columns where e4 is NaN
    for c in ["dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
              "dnase_concentration_5bp", "dnase_max_pos_rel", "dnase_max_value"]:
        if c + "_e3" in dnase.columns:
            dnase[c] = dnase[c].fillna(dnase[c + "_e3"])
            dnase = dnase.drop(columns=[c + "_e3"])
        if c + "_e2" in dnase.columns:
            dnase[c] = dnase[c].fillna(dnase[c + "_e2"])
            dnase = dnase.drop(columns=[c + "_e2"])

    # ATAC features (Exp 002 only — SCN1A n=10+10)
    atac = e002[["gene", "chrom", "pos", "ref", "variant_label",
                 "atac_magnitude", "atac_frac_5bp", "atac_frac_15bp",
                 "atac_concentration_5bp", "atac_max_pos_rel", "atac_max_value"]]

    # Splicing features (Exp 001 only — SCN1A n=10+10)
    splice = e001[["gene", "chrom", "pos", "ref", "variant_label",
                   "splice_magnitude", "splice_frac_5bp", "splice_frac_15bp",
                   "splice_concentration_5bp", "splice_max_pos_rel",
                   "splice_max_value"]]

    # Outer-merge all three on the variant identity
    feat = (dnase
            .merge(atac, on=["gene", "chrom", "pos", "ref", "variant_label"],
                   how="outer")
            .merge(splice, on=["gene", "chrom", "pos", "ref", "variant_label"],
                   how="outer"))

    feat["label_binary"] = feat["variant_label"].map(
        {"pathogenic": 1, "benign": 0}
    ).astype("Int64")

    feat = feat.sort_values(["gene", "chrom", "pos"]).reset_index(drop=True)
    return feat


# ---------------------------------------------------------------------------
# 2. Evaluation utilities
# ---------------------------------------------------------------------------

def evaluate(y_true: np.ndarray, y_score: np.ndarray) -> dict:
    """AUROC, AUPRC, top-5% precision. y_score need not be in [0,1]."""
    if len(set(y_true.tolist())) < 2:
        return {"auroc": float("nan"), "auprc": float("nan"),
                "top5_precision": float("nan")}
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)
    n = len(y_true)
    k = max(1, int(round(0.05 * n)))
    top_k_idx = np.argsort(y_score)[::-1][:k]
    top5_prec = float(np.mean(y_true[top_k_idx]))
    return {"auroc": float(auroc), "auprc": float(auprc),
            "top5_precision": top5_prec}


def cross_validate(X: pd.DataFrame, y: np.ndarray,
                   model_factory) -> pd.DataFrame:
    """Stratified 5-fold CV. Returns per-fold metrics."""
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                          random_state=RANDOM_STATE)
    rows = []
    for fold, (tr, va) in enumerate(skf.split(X, y)):
        medians = X.iloc[tr].median(numeric_only=True)
        X_tr = X.iloc[tr].fillna(medians)
        X_va = X.iloc[va].fillna(medians)
        y_tr, y_va = y[tr], y[va]
        model = model_factory()
        model.fit(X_tr, y_tr)
        # Use decision_function when available (LogReg), else predict_proba.
        if hasattr(model, "decision_function"):
            scores = model.decision_function(X_va)
        else:
            scores = model.predict_proba(X_va)[:, 1]
        m = evaluate(y_va, scores)
        rows.append({
            "fold": fold,
            "n_train": len(tr), "n_val": len(va),
            "n_pos_train": int(y_tr.sum()), "n_pos_val": int(y_va.sum()),
            **m,
        })
    return pd.DataFrame(rows)


def per_feature_logreg_cv(X: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    """Single-feature logistic regression with the same CV scheme."""
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                          random_state=RANDOM_STATE)
    rows = []
    for col in X.columns:
        fold_auprcs, fold_aurocs, fold_top5 = [], [], []
        for tr, va in skf.split(X, y):
            med = X.iloc[tr][col].median()
            x_tr = X.iloc[tr][col].fillna(med).values.reshape(-1, 1)
            x_va = X.iloc[va][col].fillna(med).values.reshape(-1, 1)
            y_tr, y_va = y[tr], y[va]
            model = LogisticRegression(class_weight="balanced", C=1.0,
                                       max_iter=1000, random_state=RANDOM_STATE)
            model.fit(x_tr, y_tr)
            scores = model.decision_function(x_va)
            fold_auprcs.append(average_precision_score(y_va, scores))
            fold_aurocs.append(roc_auc_score(y_va, scores))
            k = max(1, int(round(0.05 * len(y_va))))
            top_k = np.argsort(scores)[::-1][:k]
            fold_top5.append(float(np.mean(y_va[top_k])))
        rows.append({
            "feature": col,
            "n_folds": len(fold_auprcs),
            "auroc_mean": float(np.mean(fold_aurocs)),
            "auroc_std":  float(np.std(fold_aurocs)),
            "auprc_mean": float(np.mean(fold_auprcs)),
            "auprc_std":  float(np.std(fold_auprcs)),
            "top5_mean":  float(np.mean(fold_top5)),
            "top5_std":   float(np.std(fold_top5)),
        })
    return pd.DataFrame(rows)


def raw_score_cv(X: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    """Per-feature raw score (no model): just use the value as score."""
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                          random_state=RANDOM_STATE)
    rows = []
    for col in X.columns:
        fold_auprcs, fold_aurocs = [], []
        for tr, va in skf.split(X, y):
            med = X.iloc[tr][col].median()
            x_va = X.iloc[va][col].fillna(med).values
            y_va = y[va]
            fold_auprcs.append(average_precision_score(y_va, x_va))
            fold_aurocs.append(roc_auc_score(y_va, x_va))
        rows.append({
            "feature": col,
            "auroc_mean_raw": float(np.mean(fold_aurocs)),
            "auroc_std_raw":  float(np.std(fold_aurocs)),
            "auprc_mean_raw": float(np.mean(fold_auprcs)),
            "auprc_std_raw":  float(np.std(fold_auprcs)),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. Main
# ---------------------------------------------------------------------------

def make_lr():
    return LogisticRegression(class_weight="balanced", C=1.0, max_iter=1000,
                              random_state=RANDOM_STATE)


def main():
    feat = build_feature_table()
    feat.to_csv(OUT / "features.csv", index=False)
    print(f"[features] shape={feat.shape}")
    print(f"[features] per-gene label counts:\n{feat.groupby(['gene','variant_label']).size()}")
    print(f"[features] overall pathogenic={int(feat['label_binary'].sum())}, "
          f"benign={int((1-feat['label_binary']).sum())}, "
          f"prevalence={feat['label_binary'].mean():.3f}")

    null_counts = feat.isna().sum()
    print(f"\n[features] NaN counts (any col with NaN):")
    print(null_counts[null_counts > 0])

    # Cast pos to int (might have read as float64)
    feat["pos"] = feat["pos"].astype(int)

    DNASE_FEATS  = ["dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
                    "dnase_concentration_5bp", "dnase_max_value"]
    ATAC_FEATS   = ["atac_magnitude", "atac_frac_5bp", "atac_frac_15bp",
                    "atac_concentration_5bp", "atac_max_value"]
    SPLICE_FEATS = ["splice_magnitude", "splice_frac_5bp", "splice_frac_15bp",
                    "splice_concentration_5bp", "splice_max_value"]
    ALL_FEATS = DNASE_FEATS + ATAC_FEATS + SPLICE_FEATS

    y = feat["label_binary"].astype(int).values

    # A. Multivariate all (DNase + ATAC + Splicing) on all variants
    print(f"\n[A] multivariate-all features: n={len(y)}  pos={int(y.sum())}  "
          f"neg={int((1-y).sum())}  features={len(ALL_FEATS)}")
    cv_all = cross_validate(feat[ALL_FEATS], y, make_lr)
    cv_all.to_csv(OUT / "cv_multivariate_all.csv", index=False)
    print(cv_all.to_string(index=False))

    # B. Multivariate DNase-only
    print(f"\n[B] DNase-only multivariate: features={len(DNASE_FEATS)}")
    cv_dnase = cross_validate(feat[DNASE_FEATS], y, make_lr)
    cv_dnase.to_csv(OUT / "cv_multivariate_dnase_only.csv", index=False)
    print(cv_dnase.to_string(index=False))

    # C. Per-feature logistic regression (on the same X_all)
    print("\n[C] per-feature logistic regression:")
    pf = per_feature_logreg_cv(feat[ALL_FEATS], y)
    pf = pf.sort_values("auprc_mean", ascending=False).reset_index(drop=True)
    pf.to_csv(OUT / "per_feature_logreg_cv.csv", index=False)
    print(pf.to_string(index=False))

    # D. Per-feature raw score
    print("\n[D] per-feature raw score (no model):")
    raw = raw_score_cv(feat[ALL_FEATS], y)
    raw = raw.sort_values("auprc_mean_raw", ascending=False).reset_index(drop=True)
    raw.to_csv(OUT / "per_feature_raw_cv.csv", index=False)
    print(raw.to_string(index=False))

    # E. SPLICE_SITES_score baseline (joined to our variants)
    # SCN1A benchmark uses 'clnsig_category' (pathogenic/benign); the others
    # use 'label' (positive/negative). Map both to binary.
    # Also normalize chrom format: SCN1A file uses 'chr2', cross-disease files
    # use 2/7/X. Coerce everything to 'chrN' string for clean joining.
    def _normalize_chrom(s):
        s = str(s)
        return s if s.startswith("chr") else "chr" + s

    spl_rows = []
    for gene, path in SPLICE_CSVS.items():
        d = pd.read_csv(path)
        d["gene"] = gene
        if "label" not in d.columns:
            d["label"] = d["clnsig_category"]
        d["chrom"] = d["chrom"].map(_normalize_chrom)
        spl_rows.append(d[["gene", "chrom", "pos", "ref", "alt",
                           "SPLICE_SITES_score", "label"]])
    spl = pd.concat(spl_rows, ignore_index=True)
    spl["pos"] = spl["pos"].astype(int)
    spl["label_binary"] = spl["label"].map(
        {"positive": 1, "pathogenic": 1, "negative": 0, "benign": 0}
    ).astype(int)

    # Aggregate SPLICE_SITES_score per (gene, chrom, pos, ref, label) — mean over alts
    spl_agg = (spl.groupby(["gene", "chrom", "pos", "ref", "label_binary"],
                           as_index=False)
                  .agg({"SPLICE_SITES_score": "mean"}))

    join = feat[["gene", "chrom", "pos", "ref", "label_binary"]].merge(
        spl_agg, on=["gene", "chrom", "pos", "ref", "label_binary"], how="left"
    )
    n_matched = int(join["SPLICE_SITES_score"].notna().sum())
    print(f"\n[E] SPLICE_SITES join: matched {n_matched} / {len(join)} variants")

    spl_folds = pd.DataFrame()
    if n_matched > 10:
        m = join.dropna()
        y_m = m["label_binary"].astype(int).values
        x_m = m["SPLICE_SITES_score"].values
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                              random_state=RANDOM_STATE)
        rows = []
        for tr, va in skf.split(m, y_m):
            y_va = y_m[va]
            x_va = x_m[va]
            if len(set(y_va.tolist())) < 2:
                continue
            k = max(1, int(round(0.05 * len(y_va))))
            rows.append({
                "fold": len(rows),
                "n_val": len(va),
                "n_pos_val": int(y_va.sum()),
                "auroc": float(roc_auc_score(y_va, x_va)),
                "auprc": float(average_precision_score(y_va, x_va)),
                "top5_prec": float(np.mean(y_va[np.argsort(x_va)[::-1][:k]])),
            })
        spl_folds = pd.DataFrame(rows)
        spl_folds.to_csv(OUT / "cv_splice_sites_baseline.csv", index=False)
        print(spl_folds.to_string(index=False))

    # F. Random baseline
    rng = np.random.default_rng(RANDOM_STATE)
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                          random_state=RANDOM_STATE)
    rand_rows = []
    for tr, va in skf.split(feat[ALL_FEATS], y):
        y_va = y[va]
        scores = rng.uniform(size=len(y_va))
        k = max(1, int(round(0.05 * len(y_va))))
        rand_rows.append({
            "fold": len(rand_rows),
            "auroc": float(roc_auc_score(y_va, scores)),
            "auprc": float(average_precision_score(y_va, scores)),
            "top5_prec": float(np.mean(y_va[np.argsort(scores)[::-1][:k]])),
        })
    rand_folds = pd.DataFrame(rand_rows)
    rand_folds.to_csv(OUT / "cv_random_baseline.csv", index=False)
    print("\n[F] random baseline:")
    print(rand_folds.to_string(index=False))

    # G. Save combined per-fold results
    cv_all.assign(model="multivariate_all").to_csv(
        OUT / "cv_combined.csv", index=False, mode="w"
    )
    for tag, df in [
        ("multivariate_dnase_only", cv_dnase),
        ("per_feature_logreg", pf.assign(model="per_feature_logreg")),
        ("per_feature_raw", raw.assign(model="per_feature_raw")),
        ("splice_sites_baseline", spl_folds.assign(model="splice_sites_baseline")),
        ("random_baseline", rand_folds.assign(model="random_baseline")),
    ]:
        df.assign(model=tag).to_csv(OUT / "cv_combined.csv",
                                    index=False, mode="a", header=False)

    # ---- H. Summary JSON ----
    def mean_std(d, c):
        return f"{d[c].mean():.3f} ± {d[c].std():.3f}"

    # Prevalence (used as the random AUPRC baseline)
    prevalence = float(y.mean())

    summary = {
        "n_variants_total": int(len(feat)),
        "n_pathogenic": int(y.sum()),
        "n_benign": int((1-y).sum()),
        "prevalence": prevalence,
        "per_gene_label_counts": (
            feat.groupby(["gene", "variant_label"]).size().unstack(fill_value=0)
            .to_dict()
        ),
        "random_baseline": {
            "auroc": mean_std(rand_folds, "auroc"),
            "auprc": mean_std(rand_folds, "auprc"),
            "top5_prec": mean_std(rand_folds, "top5_prec"),
            "auroc_expectation": "0.500",
            "auprc_expectation": f"{prevalence:.3f} (= prevalence)",
        },
        "multivariate_all_dnase_atac_splice": {
            "n_features": len(ALL_FEATS),
            "auroc": mean_std(cv_all, "auroc"),
            "auprc": mean_std(cv_all, "auprc"),
            "top5_prec": mean_std(cv_all, "top5_precision"),
        },
        "multivariate_dnase_only": {
            "n_features": len(DNASE_FEATS),
            "auroc": mean_std(cv_dnase, "auroc"),
            "auprc": mean_std(cv_dnase, "auprc"),
            "top5_prec": mean_std(cv_dnase, "top5_precision"),
        },
        "per_feature_logreg": pf.set_index("feature").to_dict(orient="index"),
        "per_feature_raw": raw.set_index("feature").to_dict(orient="index"),
        "splice_sites_baseline": (
            {"auroc": mean_std(spl_folds, "auroc"),
             "auprc": mean_std(spl_folds, "auprc"),
             "top5_prec": mean_std(spl_folds, "top5_prec"),
             "n_matched": n_matched}
            if len(spl_folds) else {"matched": n_matched, "error": "too few"}
        ),
    }
    (OUT / "results.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()