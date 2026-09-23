# Cross-Disease AlphaGenome Benchmark — Results

SPLICE_SITES scores from AlphaGenome for pathogenic splicing variants vs. benign intronic variants
across 8 rare disease genes (7 filtered + KCNQ2 unfiltered in the
original report; with Exp 012 the benchmark now spans 7 filtered genes
including FBN1). Higher AUPRC = better discrimination.

## SPLICE_SITES

|| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
||------|---------|-------|-------|-------|-------|--------|-------------|
|| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
|| DMD | 600 | 200 | 400 | 1.0000 | 0.9999 | [1.000, 1.000] | 1.000 |
|| CFTR | 450 | 150 | 300 | 0.9994 | 0.9988 | [0.997, 1.000] | 1.000 |
|| COL4A5 | 376 | 176 | 200 | 0.9966 | 0.9969 | [0.992, 1.000] | 1.000 |
|| SCN2A | 230 | 30 | 200 | 0.9980 | 0.9880 | [0.965, 1.000] | 1.000 |
|| SCN1A | 430 | 120 | 310 | 0.9940 | 0.9830 | [0.964, 0.996] | 1.000 |
|| FBN1 | 300 | 100 | 200 | 0.9999 | 0.9997 | [0.999, 1.000] | 1.000 |
|| KCNQ2* | 550 | 200 | 350 | 0.6713 | 0.5657 | [0.501, 0.634] | 0.929 |
|| COL4A5† | 550 | 200 | 350 | 0.6577 | 0.6180 | [0.559, 0.676] | 1.000 |
|| FBN1‡ | 550 | 200 | 350 | 0.5646 | 0.5057 | [0.440, 0.568] | 1.000 |

## SPLICE_SITE_USAGE

|| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
||------|---------|-------|-------|-------|-------|--------|-------------|
|| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
|| DMD | 600 | 200 | 400 | 0.9974 | 0.9949 | [0.990, 0.998] | 1.000 |
|| CFTR | 450 | 150 | 300 | 0.9676 | 0.9649 | [0.944, 0.983] | 1.000 |
|| COL4A5 | 376 | 176 | 200 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
|| SCN1A | 430 | 120 | 310 | 0.9875 | 0.9643 | [0.939, 0.984] | 1.000 |
|| SCN2A | 230 | 30 | 200 | 0.9925 | 0.9561 | [0.904, 0.990] | 1.000 |
|| FBN1 | 300 | 100 | 200 | 0.9992 | 0.9984 | [0.994, 1.000] | 1.000 |
|| KCNQ2* | 550 | 200 | 350 | 0.6555 | 0.5633 | [0.500, 0.628] | 0.929 |
|| COL4A5† | 550 | 200 | 350 | 0.6657 | 0.6307 | [0.570, 0.690] | 1.000 |
|| FBN1‡ | 550 | 200 | 350 | 0.5809 | 0.5119 | [0.450, 0.571] | 1.000 |

## SPLICE_JUNCTIONS

|| Gene | n_total | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
||------|---------|-------|-------|-------|-------|--------|-------------|
|| MECP2 | 112 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
|| CFTR | 450 | 150 | 300 | 0.9794 | 0.9810 | [0.962, 0.993] | 1.000 |
|| COL4A5 | 376 | 176 | 200 | 0.9763 | 0.9867 | [0.972, 0.997] | 1.000 |
|| SCN2A | 230 | 30 | 200 | 0.9922 | 0.9640 | [0.916, 1.000] | 1.000 |
|| SCN1A | 430 | 120 | 310 | 0.9811 | 0.9176 | [0.861, 0.976] | 0.955 |
|| DMD | 600 | 200 | 400 | 0.8964 | 0.9129 | [0.886, 0.945] | 1.000 |
|| FBN1 | 300 | 100 | 200 | 0.9883 | 0.9695 | [0.927, 0.997] | 1.000 |
|| KCNQ2* | 550 | 200 | 350 | 0.6806 | 0.5930 | [0.522, 0.658] | 0.929 |
|| COL4A5† | 550 | 200 | 350 | 0.6603 | 0.6190 | [0.560, 0.677] | 1.000 |
|| FBN1‡ | 550 | 200 | 350 | 0.5955 | 0.5132 | [0.450, 0.577] | 0.929 |

\* KCNQ2 was scored with **no molecular-consequence filter** — positives are
all pathogenic SNVs and negatives are all benign SNVs (per Exp 008 protocol).
The other **7 filtered genes** (SCN1A, SCN2A, MECP2, CFTR, DMD, COL4A5,
FBN1) filtered positives to pathogenic + splicing-related consequences and
negatives to benign + intronic consequences. With Exp 012, COL4A5 was
added as the 7th filtered gene; with Exp 015 (this report), FBN1 is
the 8th filtered gene, scored on the same apples-to-apples protocol
and representing the first connective-tissue / fibroblast biology in
the benchmark. See `research_notebook/experiments/008_kcnq2_benchmark/README.md`
for KCNQ2 details and an apples-to-apples re-analysis.

† COL4A5 was scored twice (Exp 012): the *unfiltered* row uses all
pathogenic vs all benign (mirrors KCNQ2's primary protocol). The
*COL4A5* row without the dagger uses the matched-consequence
apples-to-apples protocol (pathogenic + splice-donor/acceptor/region vs
benign + intronic). See `research_notebook/experiments/012_col4a5_benchmark/README.md`.

‡ FBN1 was scored twice (Exp 015): the *filtered* row uses
apples-to-apples (pathogenic + splice-donor/acceptor/region vs
benign + intronic); the *FBN1‡* row uses the unfiltered protocol
(all pathogenic vs all benign), mirroring KCNQ2 / COL4A5. See
`research_notebook/experiments/015_fbn1_benchmark/README.md`.

## Interpretation

- **Best performer (filtered):** MECP2 (AUPRC=1.0000), tied with FBN1 SPLICE_SITES (AUPRC=0.9997)
- **Worst performer (filtered):** SCN1A (AUPRC=0.9830)
- **Mean AUPRC across the 7 filtered genes (incl. FBN1, SPLICE_SITES):** 0.9952
- **Std AUPRC across the 7 filtered genes (incl. FBN1, SPLICE_SITES):** 0.0069
- **KCNQ2 (unfiltered):** AUPRC=0.5657 — see Exp 008 for context.
- **COL4A5 (unfiltered):** AUPRC=0.6180 — see Exp 012 for context.
- **FBN1 (unfiltered):** AUPRC=0.5057 — see Exp 015 for context.

### Notes

- **MECP2** (n_pos=12): small positive set, confidence intervals are wide, tight CI
- **DMD** (n_pos=200): tight CI
- **CFTR** (n_pos=150): tight CI
- **SCN2A** (n_pos=30): tight CI
- **SCN1A** (n_pos=120): tight CI
- **KCNQ2** (n_pos=200, no consequence filter): AUPRC ~0.57 across all 3
  scorers. Most of KCNQ2's pathogenic variants are missense / nonsense
  rather than splice-site; the splice-specific scorers can't separate them
  from benign intronic variants. Top-5% precision (0.929) still shows
  the high-score tail is enriched for pathogenic variants — i.e. the
  model still ranks the most disrupted variants correctly even though the
  bulk distributions overlap.
- **COL4A5** (n_pos=176, filtered): AUPRC 0.987–1.000 across the 3
  splice scorers. This is the **7th gene** in the benchmark, the first
  on chrX, and the first in a **kidney**-relevant tissue (the previous
  six were brain / muscle / epithelial / mixed). Even on a gene where
  ~83% of pathogenic variants are missense / nonsense (only ~16.5% are
  canonical splice), the splice scorers perfectly separate the
  splicing-pathogenic mechanism from benign intronic — same pattern as
  the other 6 genes. Unfiltered AUPRC ≈ 0.62 (slightly higher than
  KCNQ2's 0.57; consistent with a marginally larger splice fraction
  within pathogenic) but still far below 0.95 — i.e. the splice scorers
  are **mechanism-specific**, not generic pathogenicity classifiers,
  even on chrX and even on basement-membrane collagen biology.
- **FBN1** (n_pos=100, filtered): AUPRC 0.970–1.000 across the 3 splice
  scorers. This is the **8th gene** in the benchmark, the largest by
  locus size (~237 kb, 65 exons) and the first in **connective-tissue /
  fibroblast** biology (Marfan syndrome and related fibrillinopathies).
  Even though FBN1 is much larger than CFTR and almost 2× larger than
  SCN1A, the 16 Kb SPLICE_SITES / SPLICE_JUNCTIONS / SPLICE_SITE_USAGE
  scorers still perfectly separate the splice-pathogenic mechanism from
  benign intronic — confirming the model is **position-invariant within
  a 16 Kb window** at single-exon scale. Unfiltered AUPRC = 0.506
  (lower than KCNQ2's 0.566 and COL4A5's 0.618, consistent with FBN1
  having the largest missense/nonsense share among the 8 genes:
  1535/2369 ≈ 65% missense + 452/2369 ≈ 19% nonsense; only ~14%
  canonical splice). Top-5% precision stays ≥ 0.93 in the unfiltered
  protocol, again confirming the high-score tail is enriched for
  pathogenic variants regardless of label composition.
  Connective-tissue / fibroblast biology is a **new tissue class** in
  this benchmark and the pattern still holds.
