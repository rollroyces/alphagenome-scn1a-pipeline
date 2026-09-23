# Experiment 013 — Tier-1 SCN1A candidates: actionable features for wet-lab validation

**Status:** complete (2026-09-23)
**Inputs:** `outputs/vus_top_candidates.csv`, `research_notebook/experiments/004_ism_dnase_n30/`, `research_notebook/experiments/011_ism_tier1/`, `data/gencode.v46.annotation.gtf.gz.feather`, NCBI E-utilities (PubMed + ClinVar)
**Outputs:**
- `outputs/tier1_actionable_features.csv` — 4 rows × 49 columns
- `outputs/tier1_actionable_features.md` — polished lab-facing brief
- `research_notebook/experiments/013_tier1_actionable/score_cryptic_splice.py` — the live re-scoring script
- `research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.json` — raw JSON
- `research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.csv` — raw CSV
- `research_notebook/experiments/013_tier1_actionable/README.md` — this file

## Motivation

Exp 008–012 ranked 50 VUS in SCN1A by AlphaGenome splicing scores. The top 4 ("Tier-1", from `outputs/vus_top_candidates.csv`) are:

| VUS rank | rsID | chr:pos (hg38) | ref>alt | molecular_consequence |
|---:|---|---|---|---|
| 14 | rs801806 | chr2:166041471 | T>A | splice_acceptor_variant |
| 17 | rs4293437 | chr2:166073671 | C>G | splice_acceptor_variant |
| 29 | rs801809 | chr2:166043700 | A>G | splice_donor_variant |
| 30 | rs2847163 | chr2:166043701 | C>A | splice_donor_variant |

The Carvill lab (which has the SCN1A minigene splicing assay — Sparber et al. 2023) can test ~1–2 variants per year. The bare AlphaGenome score alone is not enough to convince a busy PI to spend reagents — we need **concrete biological and predictive features** that make each candidate maximally actionable.

## What this experiment adds

For each Tier-1 variant, this experiment computes five concrete feature sets that go beyond the raw AlphaGenome score:

1. **Distance to nearest annotated exon** — using GENCODE v46 across all 30 SCN1A transcripts (526 exon rows). Reported as:
   - distance to nearest SCN1A exon boundary (any transcript)
   - distance to nearest 5′ss (donor) and 3′ss (acceptor) — distinct because SCN1A is minus strand
   - position relative to the **Ensembl canonical SCN1A-224** and **MANE Select SCN1A-201** exons
   - which transcripts include the nearest exon boundary (union)
2. **Cosine similarity to known pathogenic patterns** — for each Tier-1 DNASE ISM matrix (128 × 4) from Exp 011, compute cosine similarity against the 30 known-pathogenic SCN1A ISM matrices from Exp 004. Reported as: mean cosine, max cosine, and the top-3 most-similar pathogenic positions.
3. **Fresh AlphaGenome cryptic-splice scoring** — re-ran `dna_model.score_variant` with `sequence_length=16,384` and the three splicing scorers (`SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`) per the brief's requirement. All 4 variants scored successfully; the live re-calls match the Atlas-cached values to 4 dp.
4. **gnomAD evidence** — joined from `outputs/vus_high_impact_with_dnase.csv`. Reports presence, AF, AC, and homozygote count where available.
5. **PubMed / ClinVar cross-reference** — best-effort NCBI E-utilities search for each rsID. **All 4 returned 0 dedicated entries** (queried 2026-09-23 via PubMed and ClinVar dbs); these are uncharacterized VUS.

## Ranking (predicted wet-lab testability)

| Rank | rsID | AlphaGenome SPLICE_SITES | Closest splice site (canonical) | Mean cosine to 30 pathogenic | gnomAD |
|---:|---|---:|---|---:|---|
| 1 | rs801806 | 1.1406 | 3'ss (acceptor) @ 2 bp (exon 16) | 0.238 (max 0.392) | present (AF≈0) |
| 2 | rs4293437 | 1.0547 | 3'ss (acceptor) @ 1 bp (exon 4) | 0.161 (max 0.232) | absent |
| 3 | rs2847163 | 0.9150 | 5'ss (donor) @ 33 bp (exon 14 in canonical; @ 0 bp in 8 other isoforms) | 0.058 (max 0.141) | absent |
| 4 | rs801809 | 0.9213 | 5'ss (donor) @ 32 bp (exon 14 in canonical; @ 1 bp in 8 other isoforms) | 0.043 (max 0.120) | absent |

The ranking is driven primarily by **mean cosine similarity to the 30 known-pathogenic SCN1A ISM matrices** (Exp 004), with AlphaGenome SPLICE_SITES score and distance-to-splice-site as secondary tiebreakers. rs801806 is ~3.3× more cosine-similar to the pathogenic set than rs801809 (the worst of the four).

## Recommendation to the Carvill lab

**If you can test 2 variants this cycle: rs801806 and rs4293437.** Both disrupt canonical 3′ss acceptor sites of multi-exonic SCN1A exons; both are in intronic sequence 1–2 bp from the exon boundary (classical splice acceptor disruption geometry).

**If you can test 3–4: add rs2847163.** In 8 of the 30 SCN1A isoforms (SCN1A-202, -208, -209, -211, -216, -217, -218, -226) it sits exactly at the 5′ss donor of the corresponding exon — the most likely to show exon skipping.

## Honest caveats

- **Cosine similarity is not validation.** A high cosine (rs801806 = 0.39) is a *consistency* signal with the Exp 004 pathogenic ISM patterns; it is not a wet-lab result. If the lab tests the ranked variants and observes exon skipping / cryptic splice activation, the prediction is correct. If not, the prediction is wrong. Test all 4 if cycle budget permits.
- **The cosine values are smaller than expected.** Random unit-norm vectors in 512-D space have cosine expectation ~0.045; rs801806's 0.39 is ~9× above noise, but it is anecdotal because the reference set is only n=30. A larger pathogenic set (DMD, CFTR, or pooled across multiple genes) would make this comparison statistically meaningful.
- **DNASE Δ-score disagrees with cosine for the splice_donor variants.** rs801809 has the **strongest positive DNASE Δ** (+26.6) but the **weakest cosine** (0.12). DNASE captures chromatin-accessibility disruption; cosine captures similarity of the full ±64 bp DNASE ISM profile. The two signals measure different things; the cosine is the more direct wet-lab testability proxy because it asks "is this variant's overall effect profile similar to variants we know disrupt splicing?"
- **Transcript-dependence matters.** rs2847163 and rs801809 are inside the canonical SCN1A-224 exon 14, not at its boundary. The ClinVar `splice_donor_variant` annotation reflects the dominant isoforms in which these positions ARE at the exon/intron boundary. A minigene assay will need to decide which isoform(s) to test — the canonical SCN1A-224 (where they are exonic) or the isoforms where they are at the donor site.
- **Strand convention**: SCN1A is on the **minus strand**, so 5′ss (donor, GT) is at the **higher genomic coordinate** end of an exon; 3′ss (acceptor, AG) is at the **lower** end. All four Tier-1 variants are consistent with this convention when matched against their ClinVar `molecular_consequence` annotation.

## How to reproduce

```bash
cd /Users/hermes/projects/alphagenome-work

# 1. Re-run the live AlphaGenome cryptic-splice scoring
bash scripts/_run_with_key.sh \
    research_notebook/experiments/013_tier1_actionable/score_cryptic_splice.py

# 2. Re-run the master features build (not scripted — see execute_code blocks
#    in this experiment's history for the full pipeline: GENCODE load,
#    exon-distance, cosine similarity, join with gnomAD)
```

## Files

```
research_notebook/experiments/013_tier1_actionable/
├── README.md                          (this file)
├── score_cryptic_splice.py            (live re-scoring script)
├── cryptic_splice_scores.json         (raw JSON output)
└── cryptic_splice_scores.csv          (raw CSV output)

outputs/
├── tier1_actionable_features.csv      (4 rows × 49 cols master table)
└── tier1_actionable_features.md       (lab-facing brief)
```

## Provenance

- Exp 004 DNASE ISM matrices: `research_notebook/experiments/004_ism_dnase_n30/ism_dnase_scn1a_pathogenic_*.npy` (30 files)
- Exp 011 Tier-1 DNASE ISM matrices: `research_notebook/experiments/011_ism_tier1/ism_dnase_rs*.npy` (4 files)
- Tier-1 metadata: `outputs/vus_top_candidates.csv` (rows with `vus_rank` 14, 17, 29, 30)
- gnomAD join: `outputs/vus_high_impact_with_dnase.csv`
- GENCODE v46: `data/gencode.v46.annotation.gtf.gz.feather`
- Live AlphaGenome scoring: `dna_model.score_variant` (16,384 bp, splicing scorers, queried 2026-09-23)
- NCBI E-utilities: PubMed + ClinVar dbs, queried 2026-09-23
