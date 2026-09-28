# AlphaGenome Cross-Disease Splicing Benchmark — Generalization Methodology

**Version:** 1.0 (October 2026, after 10-gene cohort)
**Audience:** Anyone wanting to add an 11th (12th, …) gene to the benchmark or apply this pipeline to a new disease-gene pair.

This document captures the full protocol we used to benchmark 10 rare disease genes (SCN1A, SCN2A, MECP2, CFTR, DMD, KCNQ2, COL4A5, FBN1, NF1, LDLR) and the lessons learned along the way. Total cohort: **4,994 scoring calls, mean filtered AUPRC 0.9844 ± 0.0152** across all 3 AlphaGenome splice scorers.

## Table of contents

1. [Quick start](#1-quick-start)
2. [How to pick a gene](#2-how-to-pick-a-gene)
3. [How to extract variants from ClinVar](#3-how-to-extract-variants-from-clinvar)
4. [How to filter (pathogenic+splice vs benign+intronic)](#4-how-to-filter-pathogenicsplice-vs-benignintronic)
5. [How to score with the AlphaGenome API](#5-how-to-score-with-the-alphagenome-api)
6. [How to compute metrics](#6-how-to-compute-metrics)
7. [How to handle missense contamination](#7-how-to-handle-missense-contamination)
8. [Common pitfalls](#8-common-pitfalls)
9. [Time budget for a single gene](#9-time-budget-for-a-single-gene)
10. [Extension to new scorers / new atlases](#10-extension-to-new-scorers--new-atlases)
11. [Validation criteria for "the methodology generalizes"](#11-validation-criteria-for-the-methodology-generalizes)

---

## 1. Quick start

For a gene called `GENESYMBOL` on chromosome `chrN`, with MANE Select transcript `ENSTXXXXXXXX.X`:

```bash
# 1. Look up coordinates in the GENCODE v46 feather
python -c "
import pandas as pd
df = pd.read_feather('data/gencode.v46.annotation.gtf.gz.feather')
g = df[(df['gene_name'] == 'GENESYMBOL') & (df['Feature'] == 'gene')]
print(g.iloc[0][['Chromosome', 'Start', 'End', 'Strand', 'gene_id']])
t = df[(df['gene_name'] == 'GENESYMBOL') & (df['Feature'] == 'transcript') &
       (df['tag'].fillna('').str.contains('MANE_Select'))]
print(t.iloc[0][['transcript_id', 'transcript_name']])
"

# 2. Copy _nf1_extract.py / _nf1_score.py to _genesymbol_extract.py / _genesymbol_score.py
# 3. Replace GENE / CHROM / START / END / OUT_TSV constants in both files
# 4. Run extraction
python scripts/_genesymbol_extract.py

# 5. Inspect the stratified TSV — confirm ≥ 100 pathogenic+splice and ≥ 200 benign+intronic
python -c "
import pandas as pd
df = pd.read_csv('outputs/_genesymbol_stratified.tsv', sep='\t')
df['mc_l'] = df['mc'].fillna('').str.lower()
print('path+splice:', ((df['clnsig_category']=='pathogenic') &
                      (df['mc_l'].str.contains('splice_donor_variant') |
                       df['mc_l'].str.contains('splice_acceptor_variant') |
                       df['mc_l'].str.contains('splice_region_variant'))).sum())
print('ben+intronic:', ((df['clnsig_category']=='benign') &
                        df['mc_l'].str.contains('intron_variant')).sum())
"

# 6. Run scoring (~3 minutes for 300 variants)
bash scripts/_run_with_key.sh scripts/_genesymbol_score.py
```

End-to-end time: **3–5 hours** including iteration, mainly spent waiting for the 300 API calls (~3 min).

---

## 2. How to pick a gene

A gene is a good candidate if it meets **all four** of the following criteria. We used these when picking NF1 and LDLR as the 9th and 10th genes.

### 2.1 Distinct tissue

Each new gene should fill a **tissue gap** in the existing cohort. Our coverage as of 10 genes:

| Tissue class | Genes |
|--------------|-------|
| Brain | SCN1A, SCN2A, MECP2, KCNQ2 |
| Muscle | DMD |
| Epithelial / lung / mixed | CFTR |
| Kidney | COL4A5 |
| Connective tissue / fibroblast | FBN1 |
| Neural-crest / Schwann cell | NF1 |
| Liver / hepatocyte | LDLR |

Aim for **≥ 7 tissue classes in a 10-gene cohort**. Future additions: retina (RB1), breast / ovary (BRCA1/2), peripheral nerve / Schwann cell (NF1 — done), colon (APC), bone marrow (RUNX1).

### 2.2 Distinct inheritance pattern

Coverage by inheritance / mechanism as of 10 genes:

| Mechanism | Genes |
|-----------|-------|
| Autosomal dominant — gain-of-function / dominant-negative | SCN1A, SCN2A, FBN1 |
| Autosomal dominant — haploinsufficient tumor suppressor | NF1 |
| Autosomal dominant — haploinsufficient common-disease | LDLR |
| Autosomal recessive | CFTR (and rare LDLR) |
| X-linked | COL4A5, MECP2 |

Future additions: mitochondrial (MT-RNR1), imprinted (SNRPN/Prader-Willi).

### 2.3 Different gene size (try a giant if you have small/medium)

| Size class | Span | Genes |
|------------|------|-------|
| Tiny | < 10 kb | (none yet) |
| Small | 10–50 kb | LDLR (44 kb) |
| Medium | 50–250 kb | SCN1A, SCN2A, MECP2, KCNQ2, COL4A5, CFTR (190 kb), DMD (~2.2 Mb, but most introns) |
| Large | 250–500 kb | NF1 (287 kb), FBN1 (237 kb) |
| Giant | > 1 Mb | DMD (technically 2.2 Mb — already in benchmark) |

Each new addition should **differ from the closest existing gene by ≥ 2× in span**. NF1 and LDLR bracket FBN1 from above and below.

### 2.4 ClinVar variant availability

**Minimum:** ≥ 100 pathogenic splice-related SNVs AND ≥ 200 benign intronic SNVs in the gene locus. We cap at 100 / 200 for cost (≈ 300 API calls = ≈ $0.30 of quota).

Probe before committing:

```python
import pandas as pd
df = pd.read_csv('outputs/_genesymbol_stratified.tsv', sep='\t')
df['mc_l'] = df['mc'].fillna('').str.lower()
pos = ((df['clnsig_category']=='pathogenic') &
       (df['mc_l'].str.contains('splice_donor_variant') |
        df['mc_l'].str.contains('splice_acceptor_variant') |
        df['mc_l'].str.contains('splice_region_variant'))).sum()
neg = ((df['clnsig_category']=='benign') &
       df['mc_l'].str.contains('intron_variant')).sum()
print(f'{pos} pathogenic+splice, {neg} benign+intronic — go: {pos >= 100 and neg >= 200}')
```

If `pos < 100` but the gene is otherwise excellent, you can lift the cap on positives (e.g., SCN2A used `n_pos = 30` because only 30 pathogenic splice variants existed in ClinVar at the time). CIs will be wide; document this in the README.

---

## 3. How to extract variants from ClinVar

### 3.1 The recipe

For gene `GENE` with NCBI Gene ID `GENE_ID`, on chromosome `CHROM`, span `[START, END]`:

```python
import subprocess

def fetch_records(chrom, start, end):
    """Fetch ClinVar records via the tabix CLI; decode bytes as UTF-8."""
    out = subprocess.run(
        ["tabix", "data/clinvar_grch38.vcf.gz", f"{chrom}:{start}-{end}"],
        capture_output=True, check=True,
    )
    return out.stdout.decode("utf-8", errors="replace").splitlines()
```

Then for each line: split by tab, parse the INFO column (8th field) as `key=value;key=value`, filter to `GENEINFO` containing `GENE:GENE_ID`, restrict to SNVs (ref.len == 1, alt.len == 1), classify `CLNSIG` into 5 buckets.

### 3.2 Why subprocess + tabix (not pysam.TabixFile.fetch())

pysam 0.24's TabixFile iterator auto-decodes bytes as **ASCII** and crashes on ClinVar records that carry non-ASCII text in any INFO sub-field. The NF1 region's first record contains `Café-au-lait_macules` (a disease name with `é`) and explodes immediately with `UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3`. We hit this on every gene after FBN1 — there are now ~15 such records in ClinVar that touch our 10-gene cohort.

The `tabix` CLI subprocess returns raw bytes, which we decode as UTF-8 with `errors="replace"`. Same speed (0.03 s for a 200 kb region), no crashes. **Use this pattern for any new gene.** See `scripts/_nf1_extract.py` and `scripts/_ldlr_extract.py` for the canonical implementation.

### 3.3 Output schema

Write to `outputs/_<gene>_stratified.tsv` with columns:

```
chrom  pos  ref  alt  rsid  clnsig  clnsig_category  mc
```

Where `clnsig_category` ∈ {pathogenic, benign, uncertain, conflicting, other} and `mc` is the raw `MC=` value (Sequence Ontology-formatted like `SO:0001575|splice_donor_variant`).

### 3.4 Coordinate convention

We use the **gene-level row** from GENCODE (not the MANE Select transcript's span). The MANE Select is typically 0–1000 bp shorter at the 5' end than the gene, and we want to capture promoter / upstream variants in the extraction (even though they almost never end up in our 100 / 200 cap after filtering).

For NF1, gene span = `chr17:31,094,926–31,382,116`; MANE Select = `chr17:31,094,976–31,377,675`. We use the gene span (50 bp wider at 5'). Same for LDLR.

---

## 4. How to filter (pathogenic+splice vs benign+intronic)

The apples-to-apples (filtered) protocol is the **primary** protocol for all 10 genes. We use this for cross-gene comparison; the unfiltered protocol was a one-time sanity check that confirmed the mechanism-specific nature of splice scorers (KCNQ2 / COL4A5 / FBN1 unfiltered AUPRC ≈ 0.51–0.62; FBN1 has the lowest splice-fraction at 13.6%, NF1 has the highest at 26%).

### 4.1 The recipe

```python
df = pd.read_csv("outputs/_<gene>_stratified.tsv", sep="\t")
df["mc_l"] = df["mc"].fillna("").str.lower()
pos = df[(df["clnsig_category"] == "pathogenic")
         & (df["mc_l"].str.contains("splice_donor_variant")
            | df["mc_l"].str.contains("splice_acceptor_variant")
            | df["mc_l"].str.contains("splice_region_variant"))].copy()
neg = df[(df["clnsig_category"] == "benign")
         & (df["mc_l"].str.contains("intron_variant"))].copy()
np.random.seed(42)
pos = pos.sample(n=min(100, len(pos)), random_state=42)
neg = neg.sample(n=min(200, len(neg)), random_state=42)
```

The substring match works because ClinVar's MC field uses `SO:XXXXXXXX|name` format; matching on `name` ignores the SO prefix.

### 4.2 Why this filter

We restrict positives to **pathogenic splice-disrupting variants** so that the splice scorers have a chance to discriminate them from benign intronic controls. If we include all pathogenic (including the ~70–85% that are missense / nonsense), the splice scores for the pathogenic bulk look just like benign intronic scores — distributions overlap and AUPRC collapses to ~0.5 (the KCNQ2 / COL4A5 / FBN1 unfiltered pattern). The filtered protocol isolates the splice-pathogenic mechanism.

### 4.3 Caps: why 100 / 200

300 variants ≈ 3 minutes of API time ≈ $0.30 of quota. AUPRC CIs tighten rapidly with n_pos ≥ 100 (see KCNQ2 n=46 vs DMD n=200 CIs). If you have < 100 pathogenic+splice variants available, drop the cap and document the wider CI in the README (SCN2A has only 30 — accepted).

### 4.4 Sanity check the filter

Before scoring, verify both pools have ≥ 50 variants and a reasonable ratio (pos:neg between 1:5 and 1:1):

```python
print(f'pos={len(pos)}, neg={len(neg)}, ratio={len(pos)/max(1,len(neg)):.2f}')
```

If ratio < 1:5, the AUPRC baseline shifts and cross-gene comparison is unfair. If ratio > 1:1, you've oversampled positives.

---

## 5. How to score with the AlphaGenome API

### 5.1 The recipe

```python
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

dna_model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

variant = genome.Variant(
    chromosome="chrN", position=int(pos),
    reference_bases=str(ref), alternate_bases=str(alt),
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[name]
           for name in ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]]
scores = dna_model.score_variant(interval=interval, variant=variant,
                                 variant_scorers=scorers)
agg = float(scores[i].X.sum())  # one scalar per scorer, summing all tracks
```

### 5.2 Why 16 Kb window

The model's SEQUENCE_LENGTH_16KB is the smallest input that gives a meaningful splice prediction. Smaller windows (e.g. 8 Kb) lose too much intronic / exonic context for the splice site / usage / junction scorers. Larger windows (e.g. 100 Kb) dramatically increase compute cost without improving AUPRC in our experiments.

### 5.3 Why all 3 splice scorers

Each scorer measures a different splice-axis:
- **SPLICE_SITES** — direct perturbation of splice-site strength (highest AUPRC, lowest variance across genes).
- **SPLICE_SITE_USAGE** — change in relative splice-isoform usage (catches exonic splice enhancers / silencers; second-highest AUPRC).
- **SPLICE_JUNCTIONS** — change in junction-spanning reads / contact (catches deep-intronic effects and exon skipping; most variable across genes).

Report all 3 in the per-gene README. The headline number is the **mean AUPRC across the 3 scorers**, not the maximum (avoids cherry-picking).

### 5.4 Per-variant aggregation

We sum `ann.X` across all output tracks: `float(scores[i].X.sum())`. This is the simplest aggregate. Alternatives (max, mean, top-k tracks) give near-identical rankings in our experiments; sum is robust and cheap.

### 5.5 Smoke test before each run

Every score script does a single trivial call before starting the loop:

```python
v = genome.Variant(chromosome="chrN", position=<somewhere in gene>,
                   reference_bases="A", alternate_bases="G")
interval = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
dna_model.score_variant(interval=interval, variant=v,
                        variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"]])
```

If this fails, exit non-zero immediately — the loop will fail too, and you'd waste quota.

---

## 6. How to compute metrics

### 6.1 The recipe

```python
from sklearn.metrics import average_precision_score, roc_auc_score
import numpy as np

def bootstrap_auprc_ci(y_true, y_score, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y_true[idx])) < 2: continue
        aucs.append(average_precision_score(y_true[idx], y_score[idx]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

y_true = (df["label"] == "positive").astype(int).values
y_score = df["SPLICE_SITES_score"].values  # one of the 3 scorers
auroc = float(roc_auc_score(y_true, y_score))
auprc = float(average_precision_score(y_true, y_score))
ci_lo, ci_hi = bootstrap_auprc_ci(y_true, y_score)
n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
top_idx = np.argsort(y_score)[-n_top:]
top_5 = float(y_true[top_idx].mean())
```

### 6.2 Why these 4 metrics

- **AUROC** — standard, easy to interpret, but insensitive to class imbalance.
- **AUPRC** — primary metric. Sensitive to imbalance; baseline = positive rate (≈ 0.33 for our 100/200 cap).
- **95% bootstrap CI on AUPRC** — distinguishes "AUPRC = 0.99" from "AUPRC = 0.999" honestly. 1000 iterations, percentile method, seed=42.
- **Top-5% precision** — clinical translation metric. If a clinician looks at the top 5% of scored VUS, what fraction are truly pathogenic? ≥ 0.95 across all 10 genes (most = 1.000).

### 6.3 Per-scorer output schema

Save to `outputs/cross_disease_<gene>_filtered_metrics.csv` with columns:

```
gene  scorer  n_total  n_positive  n_negative  auroc  auprc  auprc_ci_lo  auprc_ci_hi  top_5_pct_precision
```

3 rows (one per scorer). The companion `outputs/cross_disease_<gene>_filtered_raw.csv` has all 300 variants × 3 scorers + metadata.

---

## 7. How to handle missense contamination

This is the **single most important methodological lesson** from the 10-gene cohort. Missense contamination is what makes the unfiltered protocol fail.

### 7.1 The phenomenon

When "pathogenic" includes the ~70–85% that are missense / nonsense, the splice scores for the pathogenic bulk look just like benign intronic scores:

| Gene | Splice-fraction of pathogenic | Unfiltered AUPRC | Filtered AUPRC |
|------|-------------------------------|------------------|----------------|
| KCNQ2 | ~16% | 0.566 | 0.990 |
| COL4A5 | ~16.5% | 0.618 | 0.995 |
| FBN1 | ~14% | 0.506 | 0.989 |
| NF1 | ~26% | (not run) | 0.998 |
| LDLR | ~12% | (not run) | 0.998 |

Pattern: lower splice-fraction → lower unfiltered AUPRC. Top-5% precision stays ≥ 0.93 in unfiltered — the high-score tail is robust regardless of label composition.

### 7.2 The lesson

**AlphaGenome's splice scorers are mechanism-specific, not generic pathogenicity classifiers.** They discriminate the splice-pathogenic mechanism from benign intronic controls. They cannot — and should not — discriminate coding pathogenic variants (missense / nonsense) from benign coding controls. That's a different prediction task (AlphaGenome has a separate `coding_score` for that).

### 7.3 What to do

1. Always report the **filtered (apples-to-apples) AUPRC** as the headline number for cross-disease comparison.
2. The unfiltered AUPRC is a useful **negative control** — if it's still ≥ 0.95, your gene is too easy (or your filtering is broken). For the splice-fraction range we've seen (12–30%), expect unfiltered AUPRC 0.50–0.65.
3. If you want a general-purpose pathogenicity predictor, score coding variants with `coding_score` and splice variants with the 3 splice scorers separately. Combine downstream.

### 7.4 What we do *not* do

We do not exclude missense variants from the "pathogenic" pool in the unfiltered run — that would be circular. The whole point of the unfiltered run is to expose the mechanism-specificity. The unfiltered AUPRC drop IS the lesson.

---

## 8. Common pitfalls

### 8.1 pysam decode error (already discussed in §3.2)

Use subprocess + tabix.

### 8.2 MANE Select transcript ID drift

The MANE Select transcript can change between GENCODE releases. Always re-verify against the GENCODE feather before scoring. For our cohort:

| Gene | MANE Select | Release |
|------|-------------|---------|
| SCN1A | ENST00000303395.9 | v46 |
| SCN2A | ENST00000283256.10 | v46 |
| MECP2 | ENST00000303391.11 | v46 |
| CFTR | ENST00000003084.11 | v46 |
| DMD | ENST00000357033.9 | v46 |
| KCNQ2 | ENST00000356457.6 | v46 |
| COL4A5 | ENST00000328300.11 | v46 |
| FBN1 | ENST00000316623.10 | v46 |
| NF1 | ENST00000358273.9 | v46 |
| LDLR | ENST00000558518.6 | v46 |

Always log which release you used.

### 8.3 Strand orientation

The benchmark is **strand-agnostic** — we pass variants to AlphaGenome as-is, with whatever strand they came from ClinVar. ClinVar reports variants on the reference (+) strand. Don't reverse-complement or flip them. We confirmed FBN1 (minus strand) gives the same results as the plus-strand genes, so strand orientation doesn't matter for this benchmark.

### 8.4 Coordinate off-by-one

ClinVar uses 1-based coordinates; the GENCODE feather also uses 1-based. pysam/tabix use 1-based. The AlphaGenome `genome.Variant` constructor uses 1-based. **No conversion needed.** If you see a 0-based / 1-based mismatch, that's a bug in your script.

### 8.5 Indel / multi-allelic handling

We restrict to SNVs (`ref.len == 1`, `alt.len == 1`). Multi-allelic records take only the first SNV alt and skip the rest. Indels and MNVs would need a separate score path (AlphaGenome supports them, but the 3 splice scorers are calibrated for SNVs in our experiments).

### 8.6 Empty positive pool

If `len(pos) < 30`, AUPRC and its CI become unreliable. If `len(pos) < 10`, don't run — pick a different gene or wait for ClinVar to grow.

### 8.7 Lifting the pos cap

If you have a gene with very few pathogenic splice variants (e.g. SCN2A has only 30), drop the cap to `min(100, len(pos))`. Document the wide CI explicitly in the README.

### 8.8 Lifting the neg cap

If you have < 200 benign+intronic variants, drop the cap. Don't go below 100 — the AUPRC baseline (random positive rate) starts to dominate the score.

### 8.9 Truncating `mc` substrings

`mc_l.str.contains('splice_donor_variant')` works because the MC field is `SO:0001575|splice_donor_variant` — matching the suffix. If ClinVar switches to a different format, update the substring. Always inspect `df['mc'].value_counts()` first.

---

## 9. Time budget for a single gene

| Phase | Time | Notes |
|-------|------|-------|
| Coordinate lookup | 5 min | 1 line of pandas on the feather |
| Copy + edit extract script | 5 min | Constant-replacement only |
| Run extraction | < 30 s | ~12,000 records / 0.03 s |
| Inspect stratified TSV | 5 min | Confirm pos+neg counts and ratio |
| Copy + edit score script | 5 min | Same structure as `_nf1_score.py` |
| Smoke test | 30 s | First API call — catches bad credentials |
| Score 300 variants | ~3 min | ~0.6 s/call × 300 |
| Compute metrics | < 5 s | Pure pandas + sklearn |
| Write README | 30 min | Mirror the NF1 / LDLR templates |
| Total | **3–5 hours** | Includes iteration if a step fails |

The dominant cost is the API calls (3 min) plus writing the README (30 min). Engineering time is small because the script templates are now well-tested.

If you have access to the AlphaGenome **Atlas** (pre-computed scores), extraction + scoring collapses to < 1 hour because you skip the live API entirely. Atlas has ClinVar variants pre-scored for GRCh38 — see `scripts/atlas_scn1a.py` for an example. Note: Atlas coverage is currently limited to ~5,000 variants per gene; not all 300 of our cap may be available for every gene.

---

## 10. Extension to new scorers / new atlases

The protocol generalizes to:

1. **Other variant scorers** — replace `SCORER_NAMES` in the score script with `["DNASE", "ATAC", "CAGE"]` to measure chromatin / promoter effects. The same `make_filtered_benchmark` and `compute_metrics` machinery applies unchanged.
2. **Other organisms / assemblies** — change `tabix data/clinvar_grch38.vcf.gz` to the GRCh37 / T2T / mouse ClinVar mirror. AlphaGenome supports GRCh38 natively; other assemblies need liftover.
3. **Other label sources** — replace ClinVar with COSMIC (somatic), gnomAD (population), or in-house functional data. The filter pattern (pathogenic + X vs benign + Y) generalizes.
4. **Stratification by exon / domain** — instead of one benchmark pool, stratify by which exon the variant falls in. Useful for detecting exon-specific model biases.

---

## 11. Validation criteria for "the methodology generalizes"

We declare the methodology **generalized** to a new gene if:

1. **AUPRC ≥ 0.95** on all 3 scorers (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) under the filtered protocol.
2. **Top-5% precision ≥ 0.95** under the filtered protocol.
3. **n_pos ≥ 30** (else CIs too wide to claim generalization).
4. **≥ 1 new tissue class** OR **≥ 1 new inheritance pattern** OR **≥ 2× size difference** from the closest existing gene.

All 10 genes in our cohort satisfy criteria 1–3. Criteria 4 was satisfied at every gene-addition step (NF1 added neural-crest tissue + haploinsufficiency; LDLR added liver tissue + common-disease biology).

The current 10-gene mean filtered AUPRC is **0.9844 ± 0.0152**, well above the 0.95 threshold. Any new gene that drops the cohort mean below 0.95 is a real signal that the methodology has limits — investigate before adding more genes.

---

## Appendix: File map

| File | Purpose |
|------|---------|
| `scripts/_<gene>_extract.py` | Pull SNVs from `data/clinvar_grch38.vcf.gz`, stratify by CLNSIG, write TSV |
| `scripts/_<gene>_score.py` | Score TSV variants via AlphaGenome live API, compute metrics |
| `outputs/_<gene>_stratified.tsv` | Full stratified variant list (1 header + N rows) |
| `outputs/cross_disease_<gene>_filtered_raw.csv` | 300 scored variants × 3 scorers + metadata |
| `outputs/cross_disease_<gene>_filtered_metrics.csv` | 3-row per-scorer metric summary |
| `outputs/cross_disease_summary.md` | Cross-disease rollup table |
| `research_notebook/experiments/<NNN>_<gene>/README.md` | Per-gene experimental write-up |
| `paper/preprint.md` §3.6 | Cross-disease generalization section in the paper |
| `docs/METHODOLOGY.md` | This document |

---

*Last updated: September 28, 2026 (after 10-gene cohort; Exp 016 NF1 + Exp 017 LDLR). For questions or to propose a new gene, open an issue on the project repo.*
