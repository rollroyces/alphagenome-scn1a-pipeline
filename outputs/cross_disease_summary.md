# Cross-Disease AlphaGenome Benchmark — Results

SPLICE_SITES scores from AlphaGenome for pathogenic splicing variants vs. benign intronic variants
across 5 rare disease genes. Higher AUPRC = better discrimination.

## SPLICE_SITES

| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
|------|---------|-------|-------|-------|-------|--------|-------------|
| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
| DMD | 600 | 200 | 400 | 1.0000 | 0.9999 | [1.000, 1.000] | 1.000 |
| CFTR | 450 | 150 | 300 | 0.9994 | 0.9988 | [0.997, 1.000] | 1.000 |
| SCN2A | 230 | 30 | 200 | 0.9980 | 0.9880 | [0.965, 1.000] | 1.000 |
| SCN1A | 430 | 120 | 310 | 0.9940 | 0.9830 | [0.964, 0.996] | 1.000 |

## SPLICE_SITE_USAGE

| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
|------|---------|-------|-------|-------|-------|--------|-------------|
| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
| DMD | 600 | 200 | 400 | 0.9974 | 0.9949 | [0.990, 0.998] | 1.000 |
| CFTR | 450 | 150 | 300 | 0.9676 | 0.9649 | [0.944, 0.983] | 1.000 |
| SCN1A | 430 | 120 | 310 | 0.9875 | 0.9643 | [0.939, 0.984] | 1.000 |
| SCN2A | 230 | 30 | 200 | 0.9925 | 0.9561 | [0.904, 0.990] | 1.000 |

## SPLICE_JUNCTIONS

| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
|------|---------|-------|-------|-------|-------|--------|-------------|
| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
| CFTR | 450 | 150 | 300 | 0.9794 | 0.9810 | [0.962, 0.993] | 1.000 |
| SCN2A | 230 | 30 | 200 | 0.9922 | 0.9640 | [0.916, 1.000] | 1.000 |
| SCN1A | 430 | 120 | 310 | 0.9811 | 0.9176 | [0.861, 0.976] | 0.955 |
| DMD | 600 | 200 | 400 | 0.8964 | 0.9129 | [0.886, 0.945] | 1.000 |

## Interpretation

- **Best performer:** MECP2 (AUPRC=1.0000)
- **Worst performer:** SCN1A (AUPRC=0.9830)
- **Mean AUPRC across genes:** 0.9939
- **Std AUPRC across genes:** 0.0079

### Notes

- **MECP2** (n_pos=12): small positive set, confidence intervals are wide, tight CI
- **DMD** (n_pos=200): tight CI
- **CFTR** (n_pos=150): tight CI
- **SCN2A** (n_pos=30): tight CI
- **SCN1A** (n_pos=120): tight CI
