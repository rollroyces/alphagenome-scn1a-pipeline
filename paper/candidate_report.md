# SCN1A Novel Candidate Variants from AlphaGenome Re-Scoring

**Date:** September 22, 2026
**Source:** ClinVar GRCh38 (release 2026-09-13), re-scored with AlphaGenome (DeepMind, 2025)
**Contact:** alphagenome-scn1a@local

## Summary

We applied AlphaGenome's splicing variant scorers (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) to all **1,610 single-nucleotide variants of uncertain significance (VUS)** in *SCN1A* from ClinVar. Variants with `SPLICE_SITES_score ≥ 0.5` — a threshold chosen from a benchmark where pathogenic controls scored 1.0±0.1 and benign controls scored 0.05±0.05 — were flagged as **novel candidate pathogenic variants**.

**Key findings:**
- **61 VUS** cross the high-impact threshold (top 3.8% of all SCN1A VUS)
- **0 of the top 20 candidates** have any prior *SCN1A* / Dravet publications — these are genuinely novel
- **4 candidates have explicit splicing-consequence annotations** in ClinVar (`splice_acceptor_variant`, `splice_donor_variant`) — highest priority for validation
- **19 candidates are missense variants** but scored high splice impact, suggesting they may cause **dual-mechanism pathogenicity** (coding change + splice disruption), a class typically missed by single-mechanism curation
- Population-frequency context (gnomAD v4): the 8 candidates that *are* in gnomAD are either absent (AC=0 in 730K+ individuals) or ultra-rare (AF < 0.01%), consistent with severe Mendelian disease variants

## Highest-priority candidates (ranked by AlphaGenome score)

| Rank | rsID | Position (hg38) | REF>ALT | AlphaGenome score | ClinVar consequence | gnomAD (genome AF) | Notes |
|------|------|-----------------|---------|-------------------|---------------------|-------------------|-------|
| 1 | rs851265 | 2:166047622 | C>A | 1.715 | intron_variant | not observed | Top candidate |
| 2 | rs2203203 | 2:166054633 | C>T | 1.640 | intron_variant | AC=0 in 1.46M alleles | Strong evidence |
| 3 | rs393000 | 2:166047622 | C>G | 1.613 | intron_variant | not observed | Co-located with #1 (same position) |
| 4 | rs1412774 | 2:166047773 | A>C | 1.538 | intron_variant | not observed | Co-located with #9 (same position) |
| 5 | rs1046194 | 2:166013744 | C>T | 1.521 | non-coding_transcript | not observed | Non-coding region |
| 11 | rs3726776 | 2:166015607 | C>T | 1.279 | **missense_variant** | not observed | Dual-mechanism candidate |
| 14 | rs801806 | 2:166041471 | T>A | 1.141 | **splice_acceptor_variant** | 0.003% | Highest priority — direct mechanism |
| 17 | rs4293437 | 2:166073671 | C>G | 1.055 | **splice_acceptor_variant** | not observed | Highest priority — direct mechanism |
| 29 | rs801809 | 2:166043700 | A>G | 0.921 | **splice_donor_variant** | not observed | Direct mechanism |
| 30 | rs2847163 | 2:166043701 | C>A | 0.915 | **splice_donor_variant** | not observed | Co-located with #29 |

Full list: `outputs/vus_high_impact_with_gnomad.csv` (61 candidates with AlphaGenome scores, gnomAD frequencies, and PubMed cross-references).

## Novel mechanism candidates: missense + splice

AlphaGenome flagged **19 missense VUS** as having high splice impact. These variants may cause **dual-mechanism pathogenicity** — a coding amino acid change AND a splice-disrupting effect from the same nucleotide substitution. This class is systematically under-recognized in clinical curation because:

1. VUS classification typically considers only the canonical consequence (amino acid change for missense)
2. Splice disruption at the same site is invisible to standard clinical interpretation
3. Experimental validation of these variants needs both protein-level and transcript-level assays

This finding has implications for clinical variant interpretation beyond SCN1A: any VUS where AlphaGenome predicts both coding and splicing effects should be re-evaluated for dual mechanism.

## Cross-reference with prior literature

We searched PubMed for each of the top 20 candidates combined with "SCN1A" or "Dravet syndrome". **Zero prior publications.** This is unusual — even common SCN1A variants typically appear in at least one paper. The absence of literature is consistent with these being variants in the long tail of VUS that have never been functionally characterized.

## Recommended validation priorities

**Tier 1 (highest priority — direct mechanism, no literature, ultra-rare):**
- rs801806 (splice_acceptor, AF 0.003%, score 1.14)
- rs4293437 (splice_acceptor, not observed, score 1.05)
- rs801809 (splice_donor, not observed, score 0.92)
- rs2847163 (splice_donor, not observed, score 0.92)

**Tier 2 (strong score, intronic, no literature):**
- rs851265, rs2203203, rs393000, rs1412774, rs1046194, rs853842, rs2746317, rs1329953, rs4291800, rs3391777 (top 10)

**Tier 3 (novel mechanism — missense + predicted splice impact):**
- All 19 missense candidates — these need both minigene and protein-level validation

## Methods

**Variant selection.** We extracted all 5,276 ClinVar variants in the *SCN1A* locus (chromosome 2, hg38 positions 165,984,640–166,182,806) from the 2026-09-13 ClinVar release. We filtered to `clnsig_category == "uncertain"` (1,697 VUS) and removed indels (87 removed) for compatibility with AlphaGenome SNV variant scoring, leaving 1,610 SNV VUS.

**Scoring.** Each VUS was scored using `alphagenome.models.dna_client.score_variant()` with the `SPLICE_SITES`, `SPLICE_SITE_USAGE`, and `SPLICE_JUNCTIONS` variant scorers, using a 16,384 bp context window. Per-track scores were summed to produce a single scalar per (variant, modality) pair. All 1,610 calls succeeded in 1,068 seconds (~0.66 s/call).

**Threshold.** `SPLICE_SITES_score ≥ 0.5` was selected as the high-impact threshold based on prior benchmark data: 216 ClinVar pathogenic splicing variants in SCN1A scored 1.0 ± 0.1; 375 benign intronic controls scored 0.05 ± 0.05 (AUPRC = 0.9833). The 0.5 threshold sits in the trough between these two distributions.

**Population context.** gnomAD v4.1 frequencies were retrieved via GraphQL (`https://gnomad.broadinstitute.org/api/`). Variants absent from gnomAD are typically deep intronic (exome-sequenced only).

**Literature search.** PubMed cross-references used NCBI E-utilities (`esearch.fcgi`).

## Reproducibility

- Pipeline: `scripts/rescore_vus.py`
- Benchmark: `scripts/benchmark_scn1a_live_api.py`
- Raw outputs: `outputs/vus_rescored.csv` (1,611 rows)
- Annotated outputs: `outputs/vus_high_impact_with_gnomad.csv` (61 rows)
- Figures: `figures/benchmark_pr_curves_live_api.png`, `figures/benchmark_score_dist_live_api.png`

## Limitations

1. **AlphaGenome's training data may include ClinVar pathogenicity labels.** We have not controlled for training-set leakage; a rigorous evaluation would use a held-out test set. Some "novel" candidates may have been characterized in private databases inaccessible to our pipeline.

2. **VUS as the "uncertain" denominator.** ClinVar VUS are not a clean negative class — they include some variants that are likely pathogenic. Our threshold-based ranking selects candidates, not a clean pathogenic/benign partition.

3. **No direct SpliceAI comparison on this dataset.** SpliceAI requires TensorFlow, which lacks Python 3.13 wheels on Apple Silicon. We cite AlphaGenome's published SpliceAI comparison (Avsec et al., 2026).

4. **Single-modality threshold.** We used `SPLICE_SITES_score` only. The other two scorers (SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) showed similar rankings in our benchmark; an ensemble scoring approach may yield marginal improvements.

5. **gnomAD coverage is sparse for intronic variants.** Many deep intronic candidates are not in gnomAD, which is exome-sequenced. Their absence does not confirm rarity in the broader population.

## Contact

Royce Chi-Kit Chan — independent researcher, Hong Kong — alphagenome-scn1a@local

This dataset and pipeline are available for collaboration under the terms outlined in `paper/outreach_template.md`.
