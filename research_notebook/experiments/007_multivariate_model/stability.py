"""
Stability analysis — repeat the multivariate CV across different random
seeds. If the AUROC/AUPRC are stable, that's evidence the result is not
just a CV-folding artifact.
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=ConvergenceWarning)

OUT = Path("/Users/hermes/projects/alphagenome-work/research_notebook/experiments/007_multivariate_model")
feat = pd.read_csv(OUT / "features.csv")

DNASE_FEATS  = ["dnase_magnitude", "dnase_frac_5bp", "dnase_frac_15bp",
                "dnase_concentration_5bp", "dnase_max_value"]
ATAC_FEATS   = ["atac_magnitude", "atac_frac_5bp", "atac_frac_15bp",
                "atac_concentration_5bp", "atac_max_value"]
SPLICE_FEATS = ["splice_magnitude", "splice_frac_5bp", "splice_frac_15bp",
                "splice_concentration_5bp", "splice_max_value"]
ALL_FEATS = DNASE_FEATS + ATAC_FEATS + SPLICE_FEATS

y = feat["label_binary"].astype(int).values
X = feat[ALL_FEATS]

SEEDS = [42, 0, 7, 17, 99, 2025, 314, 1337]

rows = []
for seed in SEEDS:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    fold_metrics = []
    for tr, va in skf.split(X, y):
        medians = X.iloc[tr].median(numeric_only=True)
        X_tr = X.iloc[tr].fillna(medians)
        X_va = X.iloc[va].fillna(medians)
        model = LogisticRegression(class_weight="balanced", C=1.0,
                                   max_iter=1000, random_state=seed)
        model.fit(X_tr, y[tr])
        scores = model.decision_function(X_va)
        fold_metrics.append({
            "auroc": roc_auc_score(y[va], scores),
            "auprc": average_precision_score(y[va], scores),
        })
    fm = pd.DataFrame(fold_metrics)
    rows.append({
        "seed": seed,
        "auroc_mean": fm["auroc"].mean(),
        "auroc_std":  fm["auroc"].std(),
        "auprc_mean": fm["auprc"].mean(),
        "auprc_std":  fm["auprc"].std(),
    })

stab = pd.DataFrame(rows)
stab.to_csv(OUT / "stability_seeds.csv", index=False)
print("=== Stability across seeds (multivariate, 15 features) ===")
print(stab.to_string(index=False))
print()
print(f"Across {len(SEEDS)} seeds:")
print(f"  AUROC: {stab['auroc_mean'].mean():.3f} ± {stab['auroc_mean'].std():.3f}  (of seed-means)")
print(f"  AUPRC: {stab['auprc_mean'].mean():.3f} ± {stab['auprc_mean'].std():.3f}  (of seed-means)")
print(f"  Range AUPRC: [{stab['auprc_mean'].min():.3f}, {stab['auprc_mean'].max():.3f}]")