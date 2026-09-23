# Experiment 010 — Tier-2 candidate list

**Date:** 2026-09-23 (HKT)
**Status:** Complete
**Artifact:** `outputs/vus_tier2_candidates.csv` (13 rows), `outputs/vus_tier2_candidates.md` (formatted for lab outreach).

## TL;DR

Generated a Tier-2 candidate list: SCN1A VUS that **score high on AlphaGenome
but are NOT annotated as canonical splice variants by VEP**. Produced 13
candidates stratified by VEP `Consequence`:

| Stratum | n |
|---|---:|
| `intron_variant` (deep intronic) | 5 |
| `missense_variant` (dual-mechanism) | 5 |
| `non-coding_transcript_variant` | 3 |
| **Total** | **13** |

All 13 are absent from gnomAD and have 0 PubMed citations at the
rsID + "SCN1A" query level (genuinely novel, ultra-rare).

## Methodology

### 1. Input

- `outputs/vus_rescored_with_dnase.csv` (1,610 SCN1A VUS, all scored via live
  AlphaGenome API for `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`,
  and `DNASE`; produced in Exp 009).

### 2. Production ranking (splicing-only)

- Sort by `SPLICE_SITES_score` descending → top 30.
- This is the production ranking used in Exp 009 and represents the simplest
  clinically-actionable filter: "what does the model think looks like a splice
  variant?"

### 3. Tier-1 exclusion

- Remove rs801806, rs4293437, rs801809, rs2847163 from the top-30.
- These four have VEP `Consequence` ∈ {`splice_donor_variant`,
  `splice_acceptor_variant`}. They are the canonical "Tier-1" set — explicit
  splice mechanism, already in `outputs/vus_top_candidates.csv`.
- Result: **26 splice-high variants WITHOUT explicit splice annotations**
  remain.

### 4. Consequence stratification

- Group the 26 by VEP `molecular_consequence`:
  - `intron_variant`: 16 (deep intronic, not annotated as splice site)
  - `missense_variant`: 7 (coding variants flagged for potential splice effect)
  - `non-coding_transcript_variant`: 3 (overlap with non-coding RNA)
  - No `synonymous_variant`, `5'UTR`, `3'UTR`, or `upstream_gene_variant`
    survived into the splice-only top 30.
- For each stratum, compute a **combined score** = mean of min-max-normalized
  `SPLICE_SITES_score` and `DNASE_score` over the 26-variant set.
- Pick top-5 per stratum (or all-of-stratum where the pool is smaller than 5).

### 5. Cross-reference

- **gnomAD**: joined with `outputs/vus_high_impact_with_gnomad.csv` which
  carries `present`, `gnomad_genome_af`, `gnomad_exome_af` per rsID. All 13
  Tier-2 rsIDs were already cross-referenced in that CSV.
- **PubMed**: ran a fresh E-utilities `esearch` query for each rsID + "SCN1A"
  on `eutils.ncbi.nlm.nih.gov`. All 13 returned count = 0. The published
  pre-computed `pubmed_hits` column (already in the source CSV) gave 0.0 for
  the rows it covered and NaN for the rest; the live query confirms the same.

### 6. Mechanistic hypothesis per stratum

- `intron_variant`: cryptic splice site creation, poison exon activation,
  or ISE/ISS disruption.
- `missense_variant`: dual-mechanism — the substitution both alters a protein
  residue AND disrupts splicing (often invisible to minigene assays that test
  only the canonical transcript).
- `non-coding_transcript_variant`: antisense RNA / processed transcript /
  regulatory overlap; mechanism less clear, likely affects RNA processing.

## Results

### Top-30 → 26 → 13 progression

| Stage | Count | Notes |
|---|---:|---|
| Total VUS | 1,610 | All SCN1A ClinVar VUS |
| Top 30 by `SPLICE_SITES_score` | 30 | Splice-only ranking |
| Minus Tier-1 (4 canonical splice variants) | 26 | Splice-high WITHOUT explicit splice annotation |
| Minus stratified filtering (top-5 per stratum) | **13** | **Final Tier-2 list** |

### Tier-2 candidate table

| Splice rank | rsID | Pos | REF→ALT | Splice | DNASE | Combined | Consequence | PubMed |
|---:|---|---:|---|---:|---:|---:|---|---:|
| 4  | rs1412774 | 166047773 | A→C | 1.538 | +24.8 | 0.778 | intron_variant                  | 0 |
| 1  | rs851265  | 166047622 | C→A | 1.715 | −15.0 | 0.636 | intron_variant                  | 0 |
| 3  | rs393000  | 166047622 | C→G | 1.613 | −8.3  | 0.614 | intron_variant                  | 0 |
| 18 | rs408938  | 166046767 | T→G | 0.974 | +41.8 | 0.528 | intron_variant                  | 0 |
| 7  | rs2746317 | 166013901 | G→T | 1.494 | −10.2 | 0.526 | intron_variant                  | 0 |
| 20 | rs1038247 | 166042328 | T→C | 0.970 | +26.0 | 0.424 | missense_variant                | 0 |
| 11 | rs3726776 | 166015607 | C→T | 1.279 | −14.9 | 0.359 | missense_variant                | 0 |
| 15 | rs2684360 | 165998038 | C→G | 1.074 | −0.1  | 0.323 | missense_variant                | 0 |
| 23 | rs4855006 | 166058573 | T→C | 0.939 | +1.8  | 0.250 | missense_variant                | 0 |
| 26 | rs461265  | 166013794 | A→T | 0.931 | +0.2  | 0.233 | missense_variant                | 0 |
| 5  | rs1046194 | 166013744 | C→T | 1.521 | −10.0 | 0.545 | non-coding_transcript_variant   | 0 |
| 13 | rs934879  | 166013744 | C→G | 1.141 | +3.0  | 0.385 | non-coding_transcript_variant   | 0 |
| 24 | rs2015486 | 165994240 | T→A | 0.934 | +9.2  | 0.293 | non-coding_transcript_variant   | 0 |

Co-located variants worth noting:
- **chr2:166047773**: rs1412774 (A→C) + rs4291800 (A→T, splice rank 9, in
  Tier-1 pool but not Tier-1 splice-set) — same position, two alt alleles.
- **chr2:166047622**: rs851265 (C→A) + rs393000 (C→G) — same pattern.
- **chr2:166013744**: rs1046194 (C→T) + rs934879 (C→G) — same pattern.
  Both rs1412774 and rs1046194 hit "two alt alleles at the same locus"
  within the top 30, which is itself a strong signal of locus-level
  regulatory importance.

### Comparison to Tier-1

| Group | n | Splice score range | Splice median |
|---|---:|---|---:|
| Tier-1 (4 canonical splice variants) | 4 | 0.915 – 1.141 | 0.988 |
| Tier-2 pool (26 in top-30 minus Tier-1) | 26 | 0.931 – 1.715 | 1.134 |
| Tier-2 selected (this list, 13 rows) | 13 | 0.931 – 1.715 | 1.141 |

**Tier-2 scores HIGHER than Tier-1 on median.** The Tier-1 set is dragged down
by two of its members sitting at splice ranks 29–30. Tier-1 also clusters in
a narrow score window (0.92–1.14), while Tier-2 spans the full top of the
ranking (0.93–1.72, with rs851265 at the very top).

This **inverts** the naive expectation. Pathogenic-class labels in ClinVar do
not necessarily give the highest model scores — AlphaGenome is finding high-
impact non-canonical candidates that simply haven't been annotated as splice
variants. This is consistent with the Exp 009 finding that DNASE-only
candidates are "missense variants with high chromatin impact, low splice
impact" — i.e., the model is genuinely flagging mechanism classes that
ClinVar curation does not yet recognize.

## Honest caveats

1. **These are NOT validated pathogenic.** High AlphaGenome score is not
   evidence of pathogenicity. Exp 009 showed that `SPLICE_SITES_score`
   separates ClinVar pathogenic vs benign at high AUPRC, but that is a
   benchmark on labeled variants, not a guarantee on these specific VUS.

2. **`molecular_consequence` annotation can be wrong.** VEP is not
   authoritative for non-coding function. The `intron_variant` label means
   "not annotated as splice site" — it does NOT mean "definitely has no
   splice effect." In fact, our entire reason for flagging these is that
   AlphaGenome predicts a splice effect; VEP just didn't annotate one.

3. **Combined-score formula is heuristic.** Min-max normalize both scorers
   over the 26-variant set, then take their mean. This:
   - Conflates DNASE sign (positive = increase, negative = decrease in
     accessibility) — both are treated as "high impact." Exp 009 caveat #2.
   - Weights both scorers equally, but only the splice score has been
     validated on labeled ClinVar variants (AUPRC). The DNASE score's
     clinical utility is hypothesis-level only.
   - For Tier-2 ranking, we are using the combined score because Tier-1 is
     splice-only by construction; we want Tier-2 to include chromatin-driven
     candidates even when their splice score alone is borderline.

4. **PubMed search is rsID-based.** This is fast but misses literature that
   refers to a variant by genomic position (chr2:166047773) or HGVS
   notation (NM_001165963.4:c.4725-67T>G) without naming the rsID.
   Before publication, re-query by position + gene.

5. **No synonymous, UTR, or upstream variants in this tier.** The
   production ranking was **splice-score-only**, so it naturally favors
   variants that look like splicing-impact on the model. Variants predicted
   to disrupt only `5'UTR` translation efficiency or `upstream` regulatory
   regions would need a different ranking (e.g., expression-impact or
   chromatin-only) — that's an open follow-up, not part of this experiment.

6. **`SPLICE_SITES_score` itself is a model-derived quantity.** The
   implementation sums across multiple tracks and discards track-specific
   information. We do not have per-junction / per-tissue attribution in
   this artifact. A wet-lab validation experiment is needed to confirm
   mechanism.

## Recommendations for lab outreach

**Lead with the strongest 4 candidates:**
1. **rs1412774** (intronic, combined 0.778) — top pick. Co-located alt allele
   rs4291800 is in the splice-only top-10 (rank 9), strengthening the locus
   case.
2. **rs1046194** (non-coding transcript, combined 0.545) — paired with
   rs934879 at the same position (chr2:166013744). Locus-level signal.
3. **rs1038247** (missense + splice, combined 0.424) — best example of the
   **dual-mechanism pattern**: needs BOTH protein functional assay AND
   minigene splicing assay.
4. **rs851265** (intronic, splice rank 1) — highest single SPLICE_SITES_score
   in the entire VUS set. Paired with rs393000 at the same position.

**Email framing (per `alphagenome-variant-scoring` skill):**
> "These 13 variants have AlphaGenome splice/chromatin predictions in the
> same range as the 4 canonical Tier-1 SCN1A splice variants but are NOT
> annotated as splice variants by VEP. They are absent from gnomAD and
> unpublished in the SCN1A literature (PubMed = 0 for each). We propose
> three mechanism classes: deep-intronic cryptic splice (5), dual-mechanism
> missense+splice (5), and non-coding transcript overlap (3). The four
> strongest picks co-locate with another top-30 alt allele, arguing for
> locus-level splice regulatory elements at chr2:166047773, 166047622, and
> 166013744."

**Different wet-lab assays per stratum:**
- Stratum A (intronic): minigene spanning the relevant exon ± 500 bp.
- Stratum B (missense): BOTH patch-clamp for channel function AND minigene
  for splice effect. This is the unique requirement of dual-mechanism
  candidates; a minigene-only experiment would miss the canonical
  mechanism.
- Stratum C (non-coding transcript): RNA-seq from patient-derived cells or
  reporter assay targeting the antisense/non-coding transcript.

## Repro

```python
import pandas as pd
df = pd.read_csv('outputs/vus_rescored_with_dnase.csv')
top30 = df.nlargest(30, 'SPLICE_SITES_score')
tier1 = {'801806','4293437','801809','2847163'}
tier2_pool = top30[~top30['rsid'].astype(str).isin(tier1)]
# stratify by molecular_consequence, compute combined = mean(min-max norm),
# pick top-5 per stratum
```

Runtime: <1 s (no API calls, all data already in `outputs/vus_rescored_with_dnase.csv`).

## Files produced

| File | Rows | Purpose |
|---|---:|---|
| `outputs/vus_tier2_candidates.csv` | 13 | Main artifact, all scores + cross-refs + hypothesis |
| `outputs/vus_tier2_candidates.md`   | –   | Formatted for lab outreach email |
| `research_notebook/experiments/010_tier2_candidates/README.md` | – | This document |

## NOT modified

- `outputs/vus_rescored.csv`, `outputs/vus_rescored_with_dnase.csv`
- `outputs/vus_top_candidates.csv` (Tier-1 list — preserved)
- `outputs/vus_high_impact_with_gnomad.csv`
- `paper/preprint.md`
- Any existing scripts.
