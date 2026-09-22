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

(entries added as we go)
