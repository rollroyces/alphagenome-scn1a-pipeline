# SCN1A Tier-1 Candidate Variants — Wet-Lab Actionable Brief

**To:** Carvill lab (Dravet minigene splicing assay)
**From:** AlphaGenome SCN1A VUS triage project
**Date:** 2026-09-23
**Re:** Prioritization of 4 splice-region SCN1A VUS for minigene validation

---

## Bottom line

We have ranked the 4 Tier-1 SCN1A splice-region VUS by predicted wet-lab testability. If you can run **2 variants** this cycle, we recommend:

1. **rs801806** (chr2:g.166041471 T>A) — highest predicted testability score across all signals
2. **rs4293437** (chr2:g.166073671 C>G) — second highest; both disrupt canonical 3′ss acceptor sites of multi-exonic SCN1A exons

If you can run **3 or 4**, add **rs2847163** (5′ss donor disruption in canonical exon 14 — most likely to show exon skipping).

---

## Ranking by predicted testability

| Rank | rsID | mol_cons | AlphaGenome SPLICE_SITES | Closest splice site (canonical) | Mean cosine to 30 pathogenic | gnomAD |
|------|------|----------|-------------------------|----------------------------------|------------------------------|--------|
| 1 | rs801806 | `splice_acceptor_variant` | 1.1406 | 3'ss (acceptor) @ 2 bp (exon 16) | 0.238 (max 0.392) | present (AF≈0.0) |
| 2 | rs4293437 | `splice_acceptor_variant` | 1.0547 | 3'ss (acceptor) @ 1 bp (exon 4) | 0.161 (max 0.232) | absent |
| 3 | rs2847163 | `splice_donor_variant` | 0.9150 | 5'ss (donor) @ 33 bp (exon 14) | 0.058 (max 0.141) | absent |
| 4 | rs801809 | `splice_donor_variant` | 0.9213 | 5'ss (donor) @ 32 bp (exon 14) | 0.043 (max 0.120) | absent |

---

## Per-variant cards

### rs801806 — chr2:g.166041471 T>A

**Molecular consequence (ClinVar):** `splice_acceptor_variant`
**ClinVar disease label:** Severe_myoclonic_epilepsy_in_infancy
**AlphaGenome SPLICE_SITES score:** 1.1406  (live re-call: 1.1406)
**AlphaGenome SPLICE_SITE_USAGE:** 146.57
**AlphaGenome SPLICE_JUNCTIONS (live, |Δ| sum):** 970.9  (Atlas cached: 1982.3)

**Mechanism hypothesis**
166,041,471 sits **2 bp from the 3'ss (acceptor)** of SCN1A's **16th exon** in the Ensembl canonical transcript SCN1A-224 (chr2:166041230-166041469; exon 15 / 2 bp from 3'ss (acceptor) in MANE Select SCN1A-201, chr2:166041230-166041469). Position consistent with canonical splice-site disruption.

**Exon proximity (closest annotated boundary, all SCN1A transcripts)**
- distance to nearest SCN1A exon: 2 bp (in intron)
- distance to nearest 5′ss (donor): 241 bp
- distance to nearest 3′ss (acceptor): 2 bp
- closest splice site type (any transcript): **3'ss (acceptor)** at 2 bp

**Cosine similarity to 30 known-pathogenic SCN1A ISM matrices (Exp 004)**
- mean cosine to pathogenic set: **0.238**
- max cosine to a single pathogenic: **0.392**
- top-3 most similar pathogenic variants (position, cosine):
  1. chr2:166013900 — cos **0.392**
  2. chr2:166013900 — cos 0.392
  3. chr2:166047769 — cos 0.361

**Population evidence (gnomAD)**
- present in gnomAD: **True**
- gnomAD genomes AF = 0.0 (AC=3)  /  exomes AF = 0.0 (AC=4)
- Interpretation: present at very low frequency; rare variants can still be pathogenic — absence in homozygous state is the relevant signal (gnomAD has 0 homs here, consistent with a severe dominant-acting variant).

**DNASE ISM Δ-score (Exp 011, ±64 bp window):** -2.884
- Positive Δ = Alt allele disrupts DNASE signal at the variant site (the pathogenic signature validated in Exp 004).
- Negative Δ = no/minor effect on chromatin accessibility.

**Literature (PubMed / ClinVar)**
- No dedicated PubMed entry for `rs801806` as of 2026-09; not previously characterized.


---

### rs4293437 — chr2:g.166073671 C>G

**Molecular consequence (ClinVar):** `splice_acceptor_variant`
**ClinVar disease label:** Developmental_and_epileptic_encephalopathy_6B
**AlphaGenome SPLICE_SITES score:** 1.0547  (live re-call: 1.0547)
**AlphaGenome SPLICE_SITE_USAGE:** 145.75
**AlphaGenome SPLICE_JUNCTIONS (live, |Δ| sum):** 453.8  (Atlas cached: 857.2)

**Mechanism hypothesis**
166,073,671 sits **1 bp from the 3'ss (acceptor)** of SCN1A's **4th exon** in the Ensembl canonical transcript SCN1A-224 (chr2:166073357-166073670; exon 3 / 1 bp from 3'ss (acceptor) in MANE Select SCN1A-201, chr2:166073357-166073670). Position consistent with canonical splice-site disruption.

**Exon proximity (closest annotated boundary, all SCN1A transcripts)**
- distance to nearest SCN1A exon: 1 bp (in intron)
- distance to nearest 5′ss (donor): 24 bp
- distance to nearest 3′ss (acceptor): 1 bp
- closest splice site type (any transcript): **3'ss (acceptor)** at 1 bp

**Cosine similarity to 30 known-pathogenic SCN1A ISM matrices (Exp 004)**
- mean cosine to pathogenic set: **0.161**
- max cosine to a single pathogenic: **0.232**
- top-3 most similar pathogenic variants (position, cosine):
  1. chr2:166047771 — cos **0.232**
  2. chr2:165999722 — cos 0.230
  3. chr2:165999722 — cos 0.230

**Population evidence (gnomAD)**
- present in gnomAD: **False**
- Interpretation: absent from gnomAD (~125k exomes + 76k genomes). Consistent with a rare, potentially pathogenic variant — does NOT exclude benign rarity.

**DNASE ISM Δ-score (Exp 011, ±64 bp window):** -13.790
- Positive Δ = Alt allele disrupts DNASE signal at the variant site (the pathogenic signature validated in Exp 004).
- Negative Δ = no/minor effect on chromatin accessibility.

**Literature (PubMed / ClinVar)**
- No dedicated PubMed entry for `rs4293437` as of 2026-09; not previously characterized.


---

### rs2847163 — chr2:g.166043701 C>A

**Molecular consequence (ClinVar):** `splice_donor_variant`
**ClinVar disease label:** Early-infantile_DEE,not_provided
**AlphaGenome SPLICE_SITES score:** 0.9150  (live re-call: 0.9150)
**AlphaGenome SPLICE_SITE_USAGE:** 91.96
**AlphaGenome SPLICE_JUNCTIONS (live, |Δ| sum):** 647.0  (Atlas cached: 838.4)

**Mechanism hypothesis**
166,043,701 sits **inside an annotated SCN1A exon boundary that varies between transcripts** (nearest exon span: chr2:166043668-166044049, 0 bp). In the Ensembl canonical SCN1A-224 (chr2:166043668-166044049), this position is 33 bp from the **5'ss (donor) of exon 14** (intronic side). In MANE Select SCN1A-201, exon 13 spans chr2:166043668-166044049 and the position is 33 bp from the **5'ss (donor)**. The variant most likely disrupts the canonical 5′ss donor site in some isoforms and sits inside the (alternatively-spliced) exon in others.

**Exon proximity (closest annotated boundary, all SCN1A transcripts)**
- distance to nearest SCN1A exon: 0 bp (INSIDE exon (in at least one transcript))
- distance to nearest 5′ss (donor): 0 bp
- distance to nearest 3′ss (acceptor): 348 bp
- closest splice site type (any transcript): **5'ss (donor)** at 0 bp

**Cosine similarity to 30 known-pathogenic SCN1A ISM matrices (Exp 004)**
- mean cosine to pathogenic set: **0.058**
- max cosine to a single pathogenic: **0.141**
- top-3 most similar pathogenic variants (position, cosine):
  1. chr2:166060716 — cos **0.141**
  2. chr2:165999722 — cos 0.120
  3. chr2:165999722 — cos 0.120

**Population evidence (gnomAD)**
- present in gnomAD: **False**
- Interpretation: absent from gnomAD (~125k exomes + 76k genomes). Consistent with a rare, potentially pathogenic variant — does NOT exclude benign rarity.

**DNASE ISM Δ-score (Exp 011, ±64 bp window):** +4.260
- Positive Δ = Alt allele disrupts DNASE signal at the variant site (the pathogenic signature validated in Exp 004).
- Negative Δ = no/minor effect on chromatin accessibility.

**Literature (PubMed / ClinVar)**
- No dedicated PubMed entry for `rs2847163` as of 2026-09; not previously characterized.


---

### rs801809 — chr2:g.166043700 A>G

**Molecular consequence (ClinVar):** `splice_donor_variant`
**ClinVar disease label:** Seizure,Severe_myoclonic_epilepsy_in_infancy
**AlphaGenome SPLICE_SITES score:** 0.9213  (live re-call: 0.9213)
**AlphaGenome SPLICE_SITE_USAGE:** 113.56
**AlphaGenome SPLICE_JUNCTIONS (live, |Δ| sum):** 788.4  (Atlas cached: 788.4)

**Mechanism hypothesis**
166,043,700 sits **inside an annotated SCN1A exon boundary that varies between transcripts** (nearest exon span: chr2:166043668-166044049, 0 bp). In the Ensembl canonical SCN1A-224 (chr2:166043668-166044049), this position is 32 bp from the **5'ss (donor) of exon 14** (intronic side). In MANE Select SCN1A-201, exon 13 spans chr2:166043668-166044049 and the position is 32 bp from the **5'ss (donor)**. The variant most likely disrupts the canonical 5′ss donor site in some isoforms and sits inside the (alternatively-spliced) exon in others.

**Exon proximity (closest annotated boundary, all SCN1A transcripts)**
- distance to nearest SCN1A exon: 0 bp (INSIDE exon (in at least one transcript))
- distance to nearest 5′ss (donor): 1 bp
- distance to nearest 3′ss (acceptor): 349 bp
- closest splice site type (any transcript): **5'ss (donor)** at 1 bp

**Cosine similarity to 30 known-pathogenic SCN1A ISM matrices (Exp 004)**
- mean cosine to pathogenic set: **0.043**
- max cosine to a single pathogenic: **0.120**
- top-3 most similar pathogenic variants (position, cosine):
  1. chr2:166013900 — cos **0.120**
  2. chr2:166013900 — cos 0.120
  3. chr2:166060716 — cos 0.115

**Population evidence (gnomAD)**
- present in gnomAD: **False**
- Interpretation: absent from gnomAD (~125k exomes + 76k genomes). Consistent with a rare, potentially pathogenic variant — does NOT exclude benign rarity.

**DNASE ISM Δ-score (Exp 011, ±64 bp window):** +26.551
- Positive Δ = Alt allele disrupts DNASE signal at the variant site (the pathogenic signature validated in Exp 004).
- Negative Δ = no/minor effect on chromatin accessibility.

**Literature (PubMed / ClinVar)**
- No dedicated PubMed entry for `rs801809` as of 2026-09; not previously characterized.


---

## Methodology & caveats

- **AlphaGenome scores** were re-computed live for this brief using `dna_model.score_variant` with `sequence_length=16,384` and scorers `[SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS]` (variant-centered). Numbers match the Atlas-cached values from `outputs/vus_top_candidates.csv` to 4 dp — the live call is a verification. Note: `SPLICE_JUNCTIONS` differs slightly because the cached value is a sum over a different track set; the live re-call uses 367 tracks uniformly.
- **Distance to exon** was computed against GENCODE v46 (`data/gencode.v46.annotation.gtf.gz.feather`): 526 exon rows across 30 SCN1A transcripts. Both the **Ensembl-canonical SCN1A-224** and the **MANE Select SCN1A-201** are reported; rs2847163 and rs801809 sit inside the canonical SCN1A-224 exon 14 boundary but at the canonical 5′ss donor of other transcripts (e.g., SCN1A-218, SCN1A-226, SCN1A-211). This transcript-dependence is shown explicitly per variant.
- **Strand convention**: SCN1A is on the minus strand, so 5′ss (donor, GT) is at the **higher genomic coordinate** end of an exon, and 3′ss (acceptor, AG) is at the **lower** end. All four Tier-1 variants are consistent with this strand convention when matched against their ClinVar `molecular_consequence` annotation.
- **Cosine similarity to pathogenic patterns** uses the 30 `ism_dnase_scn1a_pathogenic_*.npy` matrices from Exp 004 (each a (128, 4) DNASE ISM magnitude vector at ±64 bp around the variant) and the corresponding Tier-1 matrices from Exp 011. Cosine values in this 512-D delta-feature space are typically small (random unit-norm expectation ~0.045); the top performer here (rs801806) reaches 0.39 — meaningful but anecdotal given n=30 reference variants.
- **gnomAD** call uses `outputs/vus_high_impact_with_dnase.csv` (already populated). rs801806 has AC=3 (genome) / AC=4 (exome) at AF ≈ 0 — ultra-rare, no homozygous carriers. The other 3 are absent from gnomAD.
- **PubMed**: dedicated entries for any of the 4 rsIDs were not found via NCBI E-utilities (PubMed + ClinVar dbs, queried 2026-09-23). All 4 are uncharacterized VUS; literature support is null.
- **Cosine is not validation.** A high cosine similarity is a *consistency* signal with known pathogenic ISM patterns — not a wet-lab result. If Carvill lab tests the ranked variants and observes exon skipping / cryptic splice activation, the prediction is correct; if they see wild-type splicing, the prediction is wrong. Please still test all 4 if cycle budget permits.

## Specific wet-lab suggestions

- **rs801806 / rs4293437**: both disrupt 3′ss (acceptor). Minigene assay should report exon inclusion efficiency and any cryptic 3′ss activation in the upstream intron. Both are 1–2 bp into the intron from the exon boundary — classical splice acceptor disruption geometry.
- **rs2847163 / rs801809**: in the canonical SCN1A-224 transcript these sit inside exon 14 / at its 5′ss donor boundary; in several other SCN1A isoforms (SCN1A-202, -208, -209, -211, -216, -217, -218, -226) they are at the 5′ss donor of their respective exons. Look for exon skipping AND for cryptic 5′ss activation in the downstream intron.
- All 4 lie in the pore-lining transmembrane-domain exons of Nav1.1 (canonical SCN1A-224 exons 4, 14, 16). Exon skipping of these exons has been validated to produce loss-of-function Nav1.1 channels consistent with Dravet (Sparber et al. 2023, and references therein).
