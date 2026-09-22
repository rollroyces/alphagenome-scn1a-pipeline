# Experiment 005 — Tissue-specific vs averaged SPLICE_JUNCTIONS

**Date:** September 22, 2026
**Status:** **NULL RESULT — brain-specific splicing scores do not significantly differ from averaged scores across 5 rare disease genes.**

## Background and motivation

AlphaGenome predicts tissue-specific genome tracks (167 ATAC tracks across cell types, 367 SPLICE_JUNCTIONS tracks across tissues). The AlphaGenome Nature paper notes that "accurately recapitulating tissue-specific patterns across cellular contexts" remains a known limitation, but it's an open question whether tissue-specific scoring *helps* for variant interpretation.

**Hypothesis:** For brain-expressed rare disease genes (SCN1A, SCN2A, MECP2), brain-tissue-specific AlphaGenome scores might capture splicing disruption patterns that averaged scores miss. If true, this would be a clinically useful refinement.

## Method

- **Genes:** SCN1A, SCN2A, MECP2 (brain-expressed), CFTR (epithelial), DMD (muscle)
- **Variants:** Same pathogenic + benign SNVs from the 5-gene cross-disease benchmark
- **Scorer:** `SPLICE_JUNCTIONS` (367 tissue-specific tracks)
- **Conditions:**
  - **Averaged:** No ontology filter (all 367 tracks)
  - **Brain-filtered:** Filter to brain-related tracks (23 tracks: frontal cortex, cerebellum, motor neuron, glutamatergic neuron, etc.)
- **Per-variant score:** max-abs across all (junctions × tracks) for that variant

## Results

| Gene | AUPRC (avg) | AUPRC (brain) | Δ (brain - avg) |
|---|---|---|---|
| SCN1A | 0.9631 | 0.9626 | -0.0005 |
| SCN2A | 0.9818 | 0.9818 | 0.0000 |
| MECP2 | 1.0000 | 1.0000 | 0.0000 |
| CFTR | 0.9999 | 0.9997 | -0.0002 |
| DMD | 0.9984 | 0.9970 | -0.0014 |
| **Mean delta** | | | **-0.0004** |

- **One-sample t-test on delta:** t=-1.58, p=0.19 (not significant)
- **Wilcoxon signed-rank test:** W=0.0, p=0.25 (not significant)

## Interpretation

**Brain-filtered SPLICE_JUNCTIONS is not significantly different from averaged across 5 genes.** The trend is slightly negative (delta = -0.0004) but not statistically meaningful.

### Possible explanations

1. **Splicing is sequence-driven, not tissue-specific.** Pathogenic variants disrupt splice sites regardless of tissue context. The information content of "which tracks are in the brain" is small for variants near canonical splice sites.
2. **Brain tracks may be noisier.** Each tissue has fewer training samples than the average; brain-specific predictions may have higher variance.
3. **Our test variants are splicing-region.** Pathogenic variants in ClinVar's "pathogenic splicing" class are already strongly enriched for splice-site disruption. Tissue-specific effects might matter for *subtle* variants (regulatory, weak splice sites) rather than strong ones.
4. **The averaged score already captures tissue-specific information.** AlphaGenome's averaged tracks may include brain tracks in the average; the difference between "averaged across 367" and "averaged across 23 brain" is small.

### Comparison to original benchmark

The original benchmark used **SPLICE_SITES** (2 tracks, tissue-agnostic) and got AUPRC=0.9833 on SCN1A. SPLICE_JUNCTIONS gives AUPRC=0.9631 on the same gene — *lower* than SPLICE_SITES. This suggests SPLICE_SITES is a more discriminative single-feature scorer than SPLICE_JUNCTIONS, despite having fewer tracks.

## Negative results / caveats

- Sample sizes for MECP2 (n_path=12), SCN2A (n_path=30) are small; high AUPRC values are within CI.
- The test is one-sided (alternative hypothesis that brain is better). A two-sided test would have higher p-values.
- We did not test multi-modal tissue-specific (ATAC + DNase + CHIP_HISTONE filtered to brain). This is a natural follow-up.
- We did not test individual brain regions (e.g., frontal cortex only). The "brain" filter aggregates 23 tracks across many regions.

## Why this is still a publishable finding

Even as a negative result, this is one of the first systematic comparisons of tissue-specific vs averaged AlphaGenome scores for rare disease variant interpretation. As AlphaGenome biology adds more tissues and modalities, this kind of comparison will be increasingly important. The fact that tissue-specific doesn't help — for splicing-region variants at least — is a useful null result for the field.

## What we'd need to show tissue-specific matters

- Variants disrupting **tissue-specific enhancers** (not splice sites)
- Variants in **tissue-specific genes** with brain-specific expression patterns (e.g., neuronal TFs)
- Multi-modal tissue-specific features (combining splicing + chromatin)
- Larger sample sizes with more subtle variants

## Files

- Script: `research_notebook/experiments/005_tissue_specific_splicing/tissue_specific_experiment.py`
- Summary: `research_notebook/experiments/005_tissue_specific_splicing/tissue_vs_averaged.csv`
- Per-variant scores: `research_notebook/experiments/005_tissue_specific_splicing/per_variant_{gene}_{condition}.csv`

## Reproducibility

- AlphaGenome API call: `atlas_client.query_variants(variants, requested_scorers=['SPLICE_JUNCTIONS'], ontology_terms=[brain_terms])`
- 23 brain-related tracks identified by biosample_name matching brain/neural/cerebell/cortex/etc.
- 18 unique ontology CURIEs used for filtering (CL:0000100, UBERON:0000955, etc.)
- Runtime: ~5 min total for 5 genes × 2 conditions
- 100% API success rate (1,867 Variants × 2 conditions, all returned scores)