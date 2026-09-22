# AlphaGenome Interpretability Research Notebook

**Question:** What has AlphaGenome learned about splicing? Can we extract
biological principles that explain why some non-coding variants are pathogenic
and most are not?

**Start date:** September 22, 2026
**Phase:** Year 1 — Foundation

## Core hypothesis

AlphaGenome's splicing predictions encode an implicit grammar of splice site
strength, conservation, and context. If we can extract this grammar, we can:

1. Predict pathogenic non-coding variants without AlphaGenome's compute
2. Discover novel rules (e.g., "deep intronic variants near exons 4-7 in
   ion channel genes are higher risk because of...")
3. Translate findings into wet-lab hypotheses

## Methodology

Three parallel lines of investigation:

### Line 1: Feature attribution
Use gradient-based or attention-based attribution to identify which sequence
positions most influence AlphaGenome's SPLICE_SITES score for known
pathogenic variants. Compare to benign controls.

### Line 2: Motif analysis
For high-impact VUS, identify common sequence motifs in the ±50 bp around
the variant. Are these motifs known (splice sites, branch points, polypyrimidine
tracts) or novel?

### Line 3: Cross-gene generalization
Test whether AlphaGenome uses the same rules across SCN1A, DMD, CFTR, etc.
If yes → universal grammar. If no → gene-specific rules.

## Log

### 2026-09-22 — Experiment 001: ISM concentration hypothesis

**Tested:** Whether pathogenic splice-disrupting variants in SCN1A show concentrated
ISM effects at the variant position vs benign intronic variants showing diffuse effects.

**Result:** **Hypothesis NOT supported.** Pathogenic variants show similar spatial
distribution to benign; the difference is in total magnitude (which AUPRC captures).
Result is null but **infrastructure is validated** — the pipeline works end-to-end
and produced 20 ISM matrices saved for follow-up analysis.

**Decision point:** Try pattern-based features (motif content of ISM matrices) or
switch modalities (ATAC, DNase, CAGE) for the next experiment. Or pivot to
paper polish while the ISM direction matures.

See `research_notebook/experiments/001_ism_scn1a/README.md` for full writeup.
