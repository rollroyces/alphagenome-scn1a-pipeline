# Comparison Table — Exp 007 Multivariate Model

**Setup:** n=147 SNVs (60 pathogenic / 87 benign) aggregated from Exp 001/002/003/004 across SCN1A, DMD, CFTR.
Per-variant ISM features from DNase (n=147 rows, all genes), ATAC (n=17 SCN1A only), splicing (n=17 SCN1A only).
5-fold stratified CV, `random_state=42`, `LogisticRegression(class_weight='balanced', C=1.0)`. Decision-function scores.

## Headline result

| Model                              | n_features | AUROC (mean ± std) | AUPRC (mean ± std) | top-5% prec (mean ± std) |
|-----------------------------------|-----------:|-------------------:|-------------------:|--------------------------:|
| **Random baseline**               |  0         | 0.515 ± 0.048      | 0.488 ± 0.081      | 0.500 ± 0.500             |
| **SPLICE_SITES_score (production)** | 1       | **1.000 ± 0.000**  | **1.000 ± 0.000**  | **1.000 ± 0.000**         |
| Multivariate (DNase only)         |  5         | 0.776 ± 0.089      | 0.748 ± 0.096      | 1.000 ± 0.000             |
| **Multivariate (all 15)**         | 15         | **0.788 ± 0.103**  | **0.766 ± 0.108**  | 1.000 ± 0.000             |
| Best single feature (dnase_frac_5bp) | 1       | 0.736 ± 0.071      | 0.694 ± 0.083      | 1.000 ± 0.000             |
| Worst single feature (atac_concentration_5bp) | 1 | 0.486 ± 0.073 | 0.418 ± 0.048 | 0.500 ± 0.447 |

## Per-fold multivariate-all detail |
| Decision:  per-fold AUROC and AUPRC for the 15-feature multivariate model

| Fold | n_train | n_val | n_pos_train | n_pos_val | AUROC   | AUPRC   |
|-----:|--------:|------:|-----------:|---------:|--------:|--------:|
| 0    | 117     | 30    | 48         | 12       | 0.875   | 0.863   |
| 1    | 117     | 30    | 48         | 12       | 0.819   | 0.772   |
| 2    | 118     | 29    | 48         | 12       | 0.887   | 0.882   |
| 3    | 118     | 29    | 48         | 12       | 0.696   | 0.656   |
| 4    | 118     | 29    | 48         | 12       | 0.662   | 0.659   |

## Per-feature ranking (logistic regression, same CV)

| Rank | Feature                  | AUROC            | AUPRC            | top-5% prec      |
|----:|--------------------------|------------------|------------------|------------------|
| 1   | dnase_frac_5bp           | 0.736 ± 0.071    | **0.694 ± 0.083**| 1.000 ± 0.000    |
| 1   | dnase_concentration_5bp  | 0.736 ± 0.071    | 0.694 ± 0.083    | 1.000 ± 0.000    |
| 3   | dnase_frac_15bp          | 0.664 ± 0.071    | 0.620 ± 0.056    | 1.000 ± 0.000    |
| 4   | dnase_magnitude          | 0.689 ± 0.067    | 0.577 ± 0.075    | 0.400 ± 0.374    |
| 5   | dnase_max_value          | 0.680 ± 0.041    | 0.545 ± 0.056    | 0.100 ± 0.200    |
| 6   | splice_magnitude         | 0.611 ± 0.025    | 0.506 ± 0.035    | 1.000 ± 0.000    |
| 7   | splice_max_value         | 0.590 ± 0.045    | 0.500 ± 0.038    | 1.000 ± 0.000    |
| 8   | atac_max_value           | 0.563 ± 0.045    | 0.479 ± 0.031    | 0.800 ± 0.245    |
| 9   | splice_frac_15bp         | 0.548 ± 0.050    | 0.471 ± 0.042    | 0.800 ± 0.400    |
| 10  | atac_magnitude           | 0.546 ± 0.049    | 0.470 ± 0.043    | 0.800 ± 0.400    |
| 11  | splice_frac_5bp          | 0.510 ± 0.051    | 0.452 ± 0.027    | 0.700 ± 0.400    |
| 11  | splice_concentration_5bp | 0.510 ± 0.051    | 0.452 ± 0.027    | 0.700 ± 0.400    |
| 13  | atac_frac_15bp           | 0.471 ± 0.060    | 0.426 ± 0.031    | 0.500 ± 0.447    |
| 14  | atac_frac_5bp            | 0.497 ± 0.070    | 0.424 ± 0.045    | 0.500 ± 0.447    |
| 15  | atac_concentration_5bp   | 0.486 ± 0.073    | 0.418 ± 0.048    | 0.500 ± 0.447    |

dnase_frac_5bp and dnase_concentration_5bp are by construction the same quantity (within-field scale *within 5bp* / total) — they tie exactly. Same for splice_frac_5bp / splice_concentration_5bp.

## Stability (8 different random seeds)

| Metric | Mean of seed-means | range of seed-means |
|--------|--------------------|---------------------|
| AUROC  | 0.790 ± 0.019      | [0.749, 0.806]      |
| AUPRC  | 0.764 ± 0.021      | [0.728, 0.795]      |

## Per-feature raw score (no model) — sanity check

Per-feature raw-score AUPRC is identical to single-feature logistic-regression AUPRC within rounding, because a single-feature logistic regression with class_weight='balanced' and C=1.0 is monotonic in the input. The lift from "raw score" to "logistic regression" is therefore zero — what matters is the multivariate model.

## Interpretation

1. **Multivariate beats the best single feature by ~0.07 AUPRC** (0.766 vs 0.694). The standard deviation of the multivariate is 0.108, so the gap is *less than one std*. This is suggestive, not definitive.
2. **The multivariate lift comes mostly from DNase features.** ATAC and splicing features are at-or-near chance on this dataset (AUROC < 0.61) because they are missing for 130/147 variants (SCN1A-only) and only 17 SCN1A rows have them.
3. **SPLICE_SITES_score is the unbeatable production baseline** at 1.0/1.0/1.0 — but this is **circular**. The 147 variants were *selected* as the top/bottom by SPLICE_SITES_score from the cross-disease benchmark files. Of course SPLICE_SITES_score perfectly separates them. The "ISM-derived multivariate" vs "SPLICE_SITES_score" comparison would only be fair on a held-out test set not used for selection — and that test set is exactly what production inference does. So SPLICE_SITES=1.0 here is the *training-set ceiling*, not a true generalization number.
4. **Top-5% precision is saturated** at 1.000 for almost every model — the dataset is too small for the tail to behave stochastically (n_val ≈ 30 per fold, top-5% = 1-2 variants).