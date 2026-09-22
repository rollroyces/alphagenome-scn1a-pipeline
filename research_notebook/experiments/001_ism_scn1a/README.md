# Experiment 001 — ISM concentration hypothesis (SCN1A)

**Date:** September 22, 2026
**Status:** **PARTIAL / INCONCLUSIVE on primary hypothesis, validated as infrastructure**

## Background and motivation

In-silico mutagenesis (ISM) systematically mutates every base in a sequence and measures the effect on model predictions. For AlphaGenome's splicing scorer, this can reveal which positions matter most for a pathogenic vs benign classification.

**Primary hypothesis:** Pathogenic splice-disrupting variants' ISM profiles should show **concentrated sensitivity at the variant position itself** (positions ±1, ±2 around the variant), because AlphaGenome "sees" the splice-site disruption locally. Benign intronic variants should show **diffuse, weaker sensitivity**, with no local concentration.

If true, this concentration metric could become a *general-purpose* feature for variant prioritization, beyond AlphaGenome's raw score.

## Method

- **Variants:** Top 10 highest-scoring pathogenic splicing variants + 10 lowest-scoring benign intronic variants from the SCN1A benchmark (Methods 3.1 of preprint).
- **Window:** 64 bp on each side of the variant (128 bp total), with 16,384 bp AlphaGenome context.
- **Scorer:** SPLICE_SITES only.
- **Metric:** Fraction of total ISM magnitude within ±5 bp and ±15 bp of the variant position.

## Results

| Group | ±5bp fraction | ±15bp fraction | Total magnitude |
|-------|---------------|----------------|-----------------|
| Pathogenic (n=10) | 0.208 ± 0.252 | 0.537 ± 0.296 | 36.5 ± 19.3 |
| Benign (n=10) | 0.116 ± 0.088 | 0.287 ± 0.172 | 11.3 ± 7.5 |

- Mann-Whitney U for ±5bp: U=50.0, **p=1.000** (not significant at α=0.05).

## Interpretation

**The primary hypothesis (concentration at variant position) is NOT supported by these data.** Pathogenic and benign variants have similar spatial distributions of ISM effects; pathogenic variants are simply larger in magnitude.

**What we learned:**

1. **ISM infrastructure works.** The pipeline (AlphaGenome `score_ism_variants` → AnnData → matrix → metrics) is now reusable.
2. **The naive concentration metric is wrong.** AlphaGenome's predictions are *not* dominated by the variant's immediate neighborhood — the model uses broader context.
3. **Magnitude differs as expected** (pathogenic variants have larger total ISM responses), but this is essentially the same signal as AUPRC already measures.

**Negative results are first-class findings.** This experiment rules out a specific, plausible hypothesis and refocuses attention on magnitude/pattern-based features.

## Next directions

The negative result suggests several productive follow-ups:

1. **Pattern-based features:** Instead of spatial concentration, look at *which sequence motifs* in the ISM response drive pathogenicity. Compare motif content between pathogenic and benign ISM matrices.
2. **Larger windows:** Try 256 bp or 1024 bp windows. AlphaGenome may use longer-range context.
3. **Other modalities:** Re-run with ATAC, DNase, or CAGE scorers. Splicing may not be the best modality for "what makes a variant pathogenic" — chromatin-level features might show cleaner separation.
4. **Other genes:** Repeat on DMD, CFTR, SCN2A to test whether the magnitude-vs-pattern relationship generalizes.

## Files

- Script: `research_notebook/experiments/001_ism_scn1a/ism_experiment.py`
- Metrics: `research_notebook/experiments/001_ism_scn1a/metrics.csv`
- ISM matrices (saved per variant): `research_notebook/experiments/001_ism_scn1a/ism_*.npy`

## Reproducibility

- AlphaGenome API call pattern: `score_ism_variants(interval, ism_interval, variant_scorers=[scorer])`
- ISM window: 64 bp on each side of variant
- Runtime: ~3 seconds per variant for 128 bp × 3 alts = 384 variant scores
- 100% API call success rate
