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

### 2026-09-22 — Experiment 001: ISM concentration hypothesis (splicing)

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

### 2026-09-22 — Experiment 002: ISM concentration in chromatin modalities

**Tested:** Whether chromatin-level features (ATAC, DNase) show a concentration
difference between pathogenic and benign SCN1A variants where splicing did not.

**Result:** **PARTIAL — DNase shows a statistically significant concentration
difference (Mann-Whitney U p=0.005, n=10+10), but ATAC does not (p=0.18).**
This is the first positive result in the ISM line. **However, n=10 is small.**
The result is hypothesis-generating, not confirmed. Replication on DMD/CFTR
is needed.

**Caveats the README understates:**
- n=10 per group is fragile; minimum two-sided p-value with this n is ~0.0001
- One outlier (benign chr2:166073346) has a huge DNase max value (217 at +32 bp)
  in an intergenic region — likely an AlphaGenome artifact, not biology
- The subagent's "real effect" framing should be tempered to "hypothesis-generating
  result pending replication"

**Honest next step:** Run on DMD and CFTR (cross-disease replication). If the
  DNase pattern holds, this becomes a real finding worth documenting in a
  methods/interpretability paper.

See `research_notebook/experiments/002_ism_atac/README.md`.

### 2026-09-22 — Experiment 003: DNase concentration replication (DMD + CFTR)

**Tested:** Whether the DNase ±5bp concentration difference from Exp 002 (SCN1A,
p=0.005) replicates on DMD and CFTR.

**Result:** **NULL — does not replicate.**

| Gene | path ±5bp | benign ±5bp | p-value |
|---|---|---|---|
| SCN1A (Exp 002) | 0.107 | 0.073 | 0.005 |
| DMD (Exp 003) | 0.112 | 0.087 | 0.137 (ns) |
| CFTR (Exp 003) | 0.082 | 0.109 | 0.594 (ns) |
| Combined (Exp 003) | 0.098 | 0.098 | 0.376 (ns) |
| Fisher meta p across genes | | | 0.286 (ns) |

**Honest interpretation:** Exp 002 was almost certainly a Type-I error from small n
(n=10) on a single test. The SCN1A "DNase concentration signature" was not a real
cross-gene effect. We have no evidence for a robust DNase-based pathogenicity feature.

**What this means:**
- The ISM concentration direction is not near-term publishable
- The methods paper (AUPRC > 0.98 across 5 genes) is unaffected — it's load-bearing
  on different evidence (splicing scores, not chromatin patterns)
- We have 60 ISM matrices saved across 3 experiments — a real dataset for future
  pattern/motif analysis
- Time to pivot: combine modalities, increase n, or try tissue-specific features

**Lesson logged:** A single p=0.005 result on n=10 is not a finding — it's a
hypothesis. Replication is the only way to know. We did the replication; it's null.
That's an honest answer, not a failure.

See `research_notebook/experiments/003_ism_dnase_replication/README.md`.
