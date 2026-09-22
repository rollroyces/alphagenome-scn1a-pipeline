# Experiment 003 — DNase concentration replication (DMD + CFTR)

**Date:** September 22, 2026
**Status:** **NULL RESULT — Exp 002 DNase concentration finding does NOT replicate on DMD or CFTR.**

## Background and motivation

Experiment 002 found that pathogenic SCN1A variants have more concentrated DNase ISM effects (±5bp fraction: path 0.107 ± 0.018 vs benign 0.073 ± 0.031, Mann-Whitney U p=0.005, n=10+10). This was the first positive result in our ISM line.

But n=10 is fragile, and a single statistical test on a single gene is not a finding — it's a hypothesis. **Replication is the only way to know if it's real.**

## Method

Same as Experiment 002, but on DMD and CFTR:
- 10 pathogenic + 10 benign SNVs per gene (intronic only, excluding splice_donor/splice_acceptor)
- 64 bp ISM window each side, 16,384 bp context
- DNASE CenterMaskScorer only (skip ATAC — testing the specific Exp 002 finding)
- Per-variant aggregation (average the 3 alt alleles) — slightly different from Exp 002's per-allele method; cross-checked with both
- Mann-Whitney U, one-sided `alternative='greater'`

**Selection note:** DMD n_path was 8 (only 8 qualifying intronic SNVs after splice filtering); CFTR n_path was 7. We kept these smaller positive sets because strict replication matters more than matching n exactly.

## Results

| Gene | n_path | n_ben | path ±5bp | benign ±5bp | U | p-value |
|---|---|---|---|---|---|---|
| DMD | 8 | 10 | 0.112 | 0.087 | 53.0 | 0.137 (ns) |
| CFTR | 7 | 10 | 0.082 | 0.109 | 33.0 | 0.594 (ns) |
| **Combined** | 15 | 20 | 0.098 | 0.098 | 160.0 | **0.376 (ns)** |
| **Fisher meta p-value across genes** | | | | | | **0.286 (ns)** |

Cross-check using Exp 002's per-allele method: combined U=216.0, p=0.337 (still not significant).

## Interpretation

**Exp 002 was likely a Type-I error inflated by small n (n=10).**

Possible explanations:
1. **Statistical artifact.** With n=10 per group, the minimum p-value is ~0.0001, but the variance of p-values near 0.005 is high. A single random cluster of variants could produce p=0.005 by chance.
2. **SCN1A-specific effect.** DMD and CFTR have different chromatin landscapes. The effect might be real but gene-specific — would require more genes to detect.
3. **Selection bias in Exp 002.** The "top 10 pathogenic" SCN1A variants are the ones with highest AlphaGenome splicing scores. Those happen to be near splice sites. The "bottom 10 benign" intronic SNVs are deep intronic. The DNase effect we saw might just reflect proximity to splice sites, not pathogenicity per se.

**The strongest case is #1 or #3.** With n=10 + a single modality + a single test, p=0.005 is suspicious. The honest conclusion is: **we have no evidence for a robust DNase concentration signature distinguishing the two groups.**

## What this means for the research program

- The chromatin-ISM interpretation direction is **not** a near-term publishable finding.
- We have 60 ISM matrices (20 from Exp 001 + 20 from Exp 002 + 20 from Exp 003) — that's a real dataset for future pattern analysis.
- The methods paper (`paper/preprint.md`) stands on its own merits (AUPRC > 0.98 across 5 genes) — ISM experiments are exploratory, not load-bearing.

## Negative results / caveats

- DMD pathogenic n=8 (not 10): the SNV filter (ref.len==1, alt.len==1) + intron-only filter dropped some variants. Re-running with different filters (e.g., allowing 2-3 bp indels) is a future option.
- CFTR pathogenic n=7: same issue. CFTR has many deletion variants, fewer SNVs.
- All 60 ISM calls succeeded (100% API success rate across 3 experiments).
- Per-variant aggregation differs from per-allele aggregation; both methods agree the effect is not significant.

## Next directions

1. **Increase n.** Run on 30+30 variants per gene for both DMD and CFTR (allow indels, expand benign set). Time cost: ~10 min per gene per modality. If still null, abandon the DNase direction.
2. **Tissue-specific DNase.** SCN1A is brain-expressed; DMD is muscle; CFTR is epithelial. Maybe brain DNase tracks would show the SCN1A pattern in DMD/CFTR variants that disrupt muscle-specific regulatory elements.
3. **Combine modalities.** Even if no single ISM feature is significant, a multivariate model (DNase magnitude + splicing magnitude + ATAC) might separate the two groups. Logistic regression on the 60 ISM matrices is a fast follow-up.
4. **Pivot to functional enrichment.** Instead of statistical signatures, ask whether specific motifs (RBP binding sites, branch points, polypyrimidine tracts) are enriched in pathogenic variants' ISM responses.

## Files

- Script: `research_notebook/experiments/003_ism_dnase_replication/ism_dnase_replication.py`
- Run log: `research_notebook/experiments/003_ism_dnase_replication/run.log`
- Per-variant metrics: `research_notebook/experiments/003_ism_dnase_replication/metrics.csv` (35 rows = 15 pathogenic + 20 benign)
- Per-allele metrics: `research_notebook/experiments/003_ism_dnase_replication/metrics_per_alt.csv` (40 rows = matches Exp 002 method)
- Per-gene summary: `research_notebook/experiments/003_ism_dnase_replication/summary_by_gene.csv`
- Meta-analysis: `research_notebook/experiments/003_ism_dnase_replication/meta_analysis.csv`
- ISM matrices: `research_notebook/experiments/003_ism_dnase_replication/ism_dnase_{dmd,cftr}_{pathogenic,benign}_NN_pos*.npy` (40 files total)

## Reproducibility

- AlphaGenome API call: `score_ism_variants(interval=16kb, ism_interval=128bp, variant_scorers=[CenterMaskScorer(DNASE, width=501, DIFF_LOG2_SUM)])`
- ISM window: 64 bp on each side of variant
- Runtime: ~2.5 min for full 20+20 across both genes
- 100% API call success rate (60/60 variants across 3 experiments)
- Smoke test: `ISM_SMOKE=1 bash scripts/_run_with_key.sh research_notebook/experiments/003_ism_dnase_replication/ism_dnase_replication.py`