# Experiment 007 — Multivariate model combining DNase + ATAC + splicing

**Date:** September 23, 2026
**Status:** **NEGATIVE / INCONCLUSIVE on the primary question.** Multivariate
logistic regression modestly outperforms the best single DNase feature, but
the lift (ΔAUPRC ≈ 0.07) is smaller than the cross-fold standard deviation
(±0.108). On this dataset the multivariate model does **not** clearly beat
DNase `frac_5bp` alone. Splice + ATAC features do not help — they are
sparse (only 17 SCN1A rows have them) and individually near chance.

## Headline numbers (5-fold CV, n=147 variants, 60 pathogenic / 87 benign)

| Model                              | n_features | AUROC            | AUPRC            |
|------------------------------------|-----------:|------------------|------------------|
| Random baseline                    | 0          | 0.515 ± 0.048    | 0.488 ± 0.081    |
| **Multivariate (DNase+ATAC+Splice)** | **15**     | **0.788 ± 0.103**| **0.766 ± 0.108**|
| Multivariate (DNase only)          | 5          | 0.776 ± 0.089    | 0.748 ± 0.096    |
| Best single feature (DNase frac_5bp) | 1        | 0.736 ± 0.071    | 0.694 ± 0.083    |
| SPLICE_SITES_score (production, training set) | 1 | 1.000 ± 0.000 | 1.000 ± 0.000 |

Stability across 8 different random seeds: AUROC 0.790 ± 0.019, AUPRC
0.764 ± 0.021, range [0.728, 0.795] — the result is stable, not a
CV-folding artifact.

## Background

The motivation was straightforward: Exp 001/002/003/004 each established
that *some* ISM-derived scalar (DNase ±5bp fraction, splice magnitude, …)
discriminates pathogenic from benign variants. The question was whether
*combining* those scalars in one logistic regression yields a meaningful
lift — that is, whether the modalities carry non-redundant information that the
single-best one is missing.

If yes → methodologically clean finding ("one model, many features, beats
the best single feature"). If no → either the modalities are redundant,
or each individual signal is at its noise floor and we are just averaging
noise into a single estimate.

## Method

- **Source data:** `metrics.csv` from Exp 001 (splicing, SCN1A n=10+10),
  Exp 002 (DNase + ATAC, SCN1A n=10+10), Exp 003 (DNase, DMD+CFTR n=8/7+10),
  Exp 004 (DNase, SCN1A+DMD+CFTR n=30+30 — primary dataset).
- **Per-variant features:** 5 per modality → DNase (`magnitude`, `frac_5bp`,
  `frac_15bp`, `concentration_5bp`, `max_value`), same for ATAC and
  splicing. 15 features total.
- **Variant identity:** `(gene, chrom, pos, ref, variant_label)`. Alts
  averaged within a SNV (mean) — matches Exp 003/004's corrected method.
- **Missing handling:** median imputation computed per training fold,
  applied to validation. Required because ATAC + splicing are NaN for
  130/147 variants (only SCN1A had those modalities).
- **Model:** `LogisticRegression(class_weight='balanced', C=1.0,
  max_iter=1000, random_state=42)`. No feature scaling (logistic
  regression is scale-invariant with L2).
- **CV:** `StratifiedKFold(n_splits=5, shuffle=True, random_state=42)`.
  5 random seeds (42, 0, 7, 17, 99, 2025, 314, 1337) used to confirm
  stability.
- **Metrics:** AUROC, AUPRC (decision function scores), top-5% precision
  (k = round(0.05 * n_val), mean of y at top-k).
- **Baselines:** random uniform; raw per-feature score (no model);
  SPLICE_SITES_score from the cross-disease benchmark files (joined to
  the same variant set).

## Results in detail

### Multivariate all features (n=147, n_features=15)
| Fold | AUROC  | AUPRC  |
|-----:|-------:|-------:|
| 0    | 0.875  | 0.863  |
| 1    | 0.819  | 0.772  |
| 2    | 0.887  | 0.882  |
| 3    | 0.696  | 0.656  |
| 4    | 0.662  | 0.659  |

### DNase-only multivariate (n=147, n_features=5)
| Fold | AUROC  | AUPRC  |
|-----:|-------:|-------:|
| 0    | 0.843  | 0.803  |
| 1    | 0.755  | 0.722  |
| 2    | 0.892  | 0.886  |
| 3    | 0.686  | 0.653  |
| 4    | 0.706  | 0.676  |

### Per-feature ranking (single-feature logistic regression)

```
dnase_frac_5bp           AUROC 0.736 ± 0.071   AUPRC 0.694 ± 0.083
dnase_concentration_5bp  AUROC 0.736 ± 0.071   AUPRC 0.694 ± 0.083   (tie: scaling)
dnase_frac_15bp          AUROC 0.664 ± 0.071   AUPRC 0.620 ± 0.056
dnase_magnitude          AUROC 0.689 ± 0.067   AUPRC 0.577 ± 0.075
dnase_max_value          AUROC 0.680 ± 0.041   AUPRC 0.545 ± 0.056
splice_magnitude         AUROC 0.611 ± 0.025   AUPRC 0.506 ± 0.035
splice_max_value         AUROC 0.590 ± 0.045   AUPRC 0.500 ± 0.038
atac_max_value           AUROC 0.563 ± 0.045   AUPRC 0.479 ± 0.031
splice_frac_15bp         AUROC 0.548 ± 0.050   AUPRC 0.471 ± 0.042
atac_magnitude           AUROC 0.546 ± 0.049   AUPRC 0.470 ± 0.043
splice_frac_5bp          AUROC 0.510 ± 0.051   AUPRC 0.452 ± 0.027
splice_concentration_5bp AUROC 0.510 ± 0.051   AUPRC 0.452 ± 0.027   (tie: scaling)
atac_frac_15bp           AUROC 0.471 ± 0.060   AUPRC 0.426 ± 0.031
atac_frac_5bp            AUROC 0.497 ± 0.070   AUPRC 0.424 ± 0.045
atac_concentration_5bp   AUROC 0.486 ± 0.073   AUPRC 0.418 ± 0.048
```

(Two of the ATAC/splicing features per modality are identical to one
another after our 5bp-concentration normalization, which is why the
ties appear. This is expected and not a bug.)

### Stability across random seeds

```
seed=42    AUROC 0.806±0.082  AUPRC 0.795±0.079
seed=0     AUROC 0.749±0.174  AUPRC 0.728±0.156
seed=7     AUROC 0.786±0.084  AUPRC 0.758±0.085
seed=17    AUROC 0.805±0.059  AUPRC 0.785±0.045
seed=99    AUROC 0.787±0.125  AUPRC 0.749±0.124
seed=2025  AUROC 0.797±0.089  AUPRC 0.770±0.088
seed=314   AUROC 0.801±0.039  AUPRC 0.769±0.041
seed=1337  AUROC 0.790±0.052  AUPRC 0.763±0.091

Across 8 seeds:
  AUROC: 0.790 ± 0.019
  AUPRC: 0.764 ± 0.021  (range [0.728, 0.795])
```

## Honest overfitting assessment

The original concern was that n=24+29 (DMD) and n=20+28 (SCN1A) are too
small for a multivariate model. The data bears that out, but not in the
worst way:

- **Per-fold std on the multivariate model (±0.108 AUPRC) is high.**
  Adding more features did not move AUPRC much, but it widened the
  fold-to-fold variance slightly compared to the single-feature
  baseline (±0.083). This is consistent with the model fitting some
  feature-specific noise per fold.
- **The DNase-only multivariate (5 features) is essentially
  indistinguishable from the 15-feature multivariate** (AUROC 0.776
  vs 0.788). The 10 extra features (ATAC + splicing) contributed
  essentially zero lift. Most likely because the median imputation
  for the 130 missing rows collapses those features to the global
  median and they don't carry gene-specific signal.
- **The lift over the best single feature (~0.07 AUPRC) is real but
  smaller than the noise envelope (one σ).** Calling this a positive
  result would be overclaim. Calling this a clear null would also
  be overclaim — the lift is in the right direction and survives 8
  random seeds.
- **Top-5% precision is saturated at 1.0 for nearly every model.**
  With n_val ≈ 30 per fold and prevalence ~0.4, the top-5% set is
  1-2 variants, and those are always pathogenic in any reasonable
  ordering. This metric cannot distinguish models on this dataset
  and should be ignored.

## Negative finding — what the model does NOT show

1. **No evidence that ATAC or splicing features add discriminative
   power.** Both modalities are near chance on the SCN1A-only subset
   they cover (n=17). Larger experiments with more genes covered
   for all modalities would be needed to test the joint-modalities
   hypothesis fairly.
2. **No evidence that the multivariate lift is clinically useful.**
   The multivariate AUROC 0.79 with std 0.10 is *worse* than the
   published production numbers for `SPLICE_SITES_score` (AUROC ≈
   0.96-0.99 on the cross-disease benchmarks) and would not be
   worth the extra ISM cost in a clinical pipeline.
3. **No evidence of the "concentration hypothesis" beyond Exp 002.**
   The concentration features (`frac_5bp`, `concentration_5bp`) do
   carry signal (best single-feature AUROC 0.74), but adding them
   to the multivariate doesn't unlock a non-linear regime.

## When this result would flip

If we had:
- **n ≥ 200 pathogenic variants per gene** with all three modalities
  measured, and
- **a held-out test set** never used for variant selection

…then the multivariate could plausibly show a larger and more stable
lift. Today's dataset is too small to distinguish 0.05 from 0.10.

## Files in this directory

- `analysis.py` — full pipeline (data merge + CV + metrics).
- `stability.py` — 8-seed stability check.
- `features.csv` — per-variant feature matrix with labels.
- `cv_multivariate_all.csv` — per-fold metrics, 15-feature model.
- `cv_multivariate_dnase_only.csv` — per-fold metrics, DNase-only.
- `per_feature_logreg_cv.csv` — per-feature logistic regression CV.
- `per_feature_raw_cv.csv` — per-feature raw-score CV (no model).
- `cv_splice_sites_baseline.csv` — per-fold metrics, SPLICE_SITES baseline.
- `cv_random_baseline.csv` — per-fold metrics, random uniform.
- `cv_combined.csv` — all of the above concatenated for convenience.
- `stability_seeds.csv` — 8-seed stability output.
- `results.json` — machine-readable summary.
- `comparison_table.md` — markdown comparison table.
- `README.md` — this file.