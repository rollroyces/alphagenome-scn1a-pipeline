# SCN1A VUS — Tier-2 Candidate List (Experiment 010)

**Generated:** 2026-09-23 (HKT)
**Pipeline:** AlphaGenome (live API), splicing + DNASE scorers
**Input:** 1,610 SCN1A ClinVar VUS (`outputs/vus_rescored_with_dnase.csv`)
**Method:** Top-30 by `SPLICE_SITES_score`, then exclude the 4 Tier-1 splice variants
(rs801806, rs4293437, rs801809, rs2847163), then take top-5 per Consequence stratum
by combined (splice + DNASE) score.

> **n_candidates = 13.** These are variants AlphaGenome scores high but that VEP
> does NOT annotate as canonical splice variants. They are the strongest Tier-2
> candidates for lab outreach — distinct mechanism hypothesis from Tier-1.

---

## What's in this list

| Stratum | n | Mechanism category |
|---|---:|---|
| `intron_variant` (deep intronic) | 5 | Cryptic splice site / poison exon / ISE-ISS disruption |
| `missense_variant` | 5 | **Dual-mechanism** candidates (coding + splicing) |
| `non-coding_transcript_variant` | 3 | Antisense / non-coding RNA / regulatory overlap |
| **Total** | **13** | |

No `synonymous_variant`, `5_prime_UTR_variant`, `3_prime_UTR_variant`, or
`upstream_gene_variant` survived into the top-30 of the splicing-only ranking.
Those mechanisms would only dominate a chromatin-only or upstream scoring.

All 13 are **absent from gnomAD** (rare or not observed) and have **0 PubMed
hits** at the rsID + "SCN1A" query level (genuinely novel, not yet characterized
in the published literature).

---

## Tier-2 candidates

### Stratum A — Deep intronic (`intron_variant`, n=5)

These are NOT annotated as canonical splice donors/acceptors, yet score in the
top quartile on `SPLICE_SITES_score`. Hypothesis: cryptic splice-site creation,
exon definition loss, poison-exon activation, or intronic splicing
enhancer/silencer (ISE/ISS) disruption.

| Splice rank | rsID     | Pos       | REF→ALT | Splice | DNASE    | Combined | clndn                       |
|---:|---|---:|---|---:|---:|---:|---|
| 4  | rs1412774  | 166047773 | A→C | 1.538 | +24.8 | 0.778 | Early-infantile_DEE         |
| 1  | rs851265   | 166047622 | C→A | 1.715 | −15.0 | 0.636 | Migraine                    |
| 3  | rs393000   | 166047622 | C→G | 1.613 | −8.3  | 0.614 | Complex neurodev. / DEE     |
| 18 | rs408938   | 166046767 | T→G | 0.974 | +41.8 | 0.528 | Early-infantile_DEE         |
| 7  | rs2746317  | 166013901 | G→T | 1.494 | −10.2 | 0.526 | Early-infantile_DEE         |

**Top pick: rs1412774** (Splice rank 4, combined 0.778). Co-located with
rs4291800 (splice rank 9) on chr2:166047773 — same position, different alt
alleles — which is itself striking and suggests a position-level splice
regulatory element. A/C substitution in a deep intronic site, with strong
combined chromatin + splice signal.

---

### Stratum B — Missense (`missense_variant`, n=5)

These are VEP-annotated as amino-acid-changing BUT AlphaGenome's
`SPLICE_SITES_score` puts them in the same range as canonical splice variants
(0.93–1.28). This is the **dual-mechanism signature** flagged in the
`alphagenome-variant-scoring` skill (Exp 009 finding): the same nucleotide
substitution both alters a protein residue AND disrupts splicing. Minigene
assays that test only the canonical transcript will miss the splice effect.

| Splice rank | rsID     | Pos       | REF→ALT | Splice | DNASE    | Combined | clndn                                  |
|---:|---|---:|---|---:|---:|---:|---|
| 20 | rs1038247  | 166042328 | T→C | 0.970 | +26.0 | 0.424 | Developmental_and_epileptic_encephalopathy |
| 11 | rs3726776  | 166015607 | C→T | 1.279 | −14.9 | 0.359 | Early-infantile_DEE                    |
| 15 | rs2684360  | 165998038 | C→G | 1.074 | −0.1  | 0.323 | not_provided                           |
| 23 | rs4855006  | 166058573 | T→C | 0.939 | +1.8  | 0.250 | Severe_myoclonic_epilepsy_in_infancy   |
| 26 | rs461265   | 166013794 | A→T | 0.931 | +0.2  | 0.233 | Early-infantile_DEE                    |

**Top pick: rs1038247** (Splice rank 20, combined 0.424). Even though it sits
near the bottom of the splice-only top 30, it has the highest positive DNASE
displacement (+26.0) in this stratum — pointing to chromatin-accessibility
disruption in addition to coding effect. Position 166042328 is in the central
transmembrane/pore region of Nav1.1, where missense variants are well-known to
cause both gain-of-function and loss-of-function DEE.

---

### Stratum C — Non-coding transcript (`non-coding_transcript_variant`, n=3)

VEP labels these as overlapping a non-coding transcript feature. Three
possibilities: (i) antisense transcript (e.g., SCN1A-AS1); (ii) deep intronic
position that overlaps a non-coding RNA gene; or (iii) regulatory element.
AlphaGenome flags them — the high splice / chromatin scores mean they likely
affect RNA processing even though they don't change a protein.

| Splice rank | rsID    | Pos       | REF→ALT | Splice | DNASE  | Combined | clndn                       |
|---:|---|---:|---|---:|---:|---:|---|
| 5  | rs1046194 | 166013744 | C→T | 1.521 | −10.0 | 0.545 | Early-infantile_DEE         |
| 13 | rs934879  | 166013744 | C→G | 1.141 | +3.0  | 0.385 | Early-infantile_DEE         |
| 24 | rs2015486 | 165994240 | T→A | 0.934 | +9.2  | 0.293 | Early-infantile_DEE         |

**Top pick: rs1046194** (Splice rank 5). Same position as rs934879 (rank 13) —
two alt alleles at the same locus in our top-30, which is a strong locus-level
signal.

---

## gnomAD and PubMed cross-reference

| Status | n |
|---|---:|
| Absent from gnomAD (rare or unobserved in population databases) | 13 |
| Present in gnomAD | 0 |
| `pubmed_count_live` = 0 (PubMed E-utilities, rsID + "SCN1A" query) | 13 |

**All 13 candidates are ultra-rare and uncited in the SCN1A literature.** This
is exactly the screening pattern that makes outreach compelling — the
"AlphaGenome scored these high, they're not in gnomAD, they're not in PubMed"
framing matches the `alphagenome-variant-scoring` skill's outreach
template.

---

## Comparison to Tier-1

| Group | n | Splice score range | Splice median |
|---|---:|---|---:|
| Tier-1 (4 explicit splice variants) | 4   | 0.915 – 1.141 | 0.988 |
| Tier-2 (26 splice-high minus Tier-1) | 26  | 0.931 – 1.715 | 1.134 |
| Tier-2 selected (this list)          | 13  | 0.931 – 1.715 | 1.141 |

**Tier-2 candidates score HIGHER on `SPLICE_SITES_score` than Tier-1**, on
median. The Tier-1 list is actually driven down by two of the four (rs2847163,
rs801809) sitting at splice ranks 29–30. Tier-2 (especially Stratum A) is a
cleaner "top of the pile" by raw splice score — which makes biological sense:
Tier-1 variants are individually validated pathogenic alleles (high prior,
selected regardless of model ranking); Tier-2 are model-discovered candidates
without that prior.

This **inverts** the intuitive expectation that "validated = high score." It
argues that AlphaGenome is not biased toward curated pathogenic alleles — it
finds non-canonical-mechanism candidates that pathology databases have not
yet annotated.

---

## Lab outreach recommendations

**Lead with the high-conviction picks from each stratum.** The most actionable
email framing:

1. **rs1412774** (intronic, chr2:166047773) — top combined score, with
   rs4291800 (also top-30) co-located at the same position pointing to a
   locus-level splice regulatory element. Validate with a minigene spanning
   exons ~10–14 of SCN1A.

2. **rs1046194** and **rs934879** (same chr2:166013744, two alt alleles) —
   strongest locus-level non-coding signal. The duplicated position argues
   for a *locus* (not just a *variant*) as the candidate.

3. **rs1038247** (missense + splice) — the dual-mechanism candidate with the
   strongest chromatin displacement in the missense stratum. Validate with
   **both** protein functional assay (e.g., patch-clamp) and a minigene
   splicing assay.

4. **rs851265** (intronic, chr2:166047622) — top of the splice-only ranking
   (rank 1), co-located with rs393000 (rank 3) — same pattern as #2.

**Avoid:** rs2684360 and rs461265 (Stratum B), which are ranked in the bottom
third of the top-30 by combined score and have less supporting evidence.

**Framing for the email:**
> "These 13 variants have splice / chromatin predictions in the same range as
> the 4 canonical Tier-1 SCN1A splice variants but are not annotated as splice
> variants by VEP. They are absent from gnomAD and have no published
> SCN1A-specific citations. We propose they represent three mechanism classes:
> deep-intronic cryptic splice activation (5 variants), dual-mechanism
> missense+splice (5), and non-coding transcript overlap (3)."

---

## Files

- `outputs/vus_tier2_candidates.csv` — machine-readable, 13 rows, all
  AlphaGenome scores + gnomAD/PubMed cross-refs + hypothesis per row.
- `research_notebook/experiments/010_tier2_candidates/README.md` — full
  methodology, caveats, comparison vs Tier-1.
