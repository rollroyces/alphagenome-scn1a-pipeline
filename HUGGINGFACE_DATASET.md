# Citable datasets and DOIs

This file tracks external citable artifacts produced from this project.

## Hugging Face Dataset

**Dataset:** RROL/scn1a-vus-alphagenome  
**URL:** https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome  
**License:** CC-BY-4.0  
**Content:** 1,610 SCN1A VUS with AlphaGenome splicing + DNASE scores, gnomAD v4.1 AF/AC/AN, PubMed hits, tier (1/2/3)  
**Citation:** Royce Lam. (2026). SCN1A VUS Re-Scoring with AlphaGenome (Revision 8db4888). Hugging Face. https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome  

### Schema (30 columns)

- `chrom, pos, ref, alt, rsid, molecular_consequence, vus_rank, tier, clnsig_category, DNASE_score, SPLICE_SITES_score, SPLICE_SITE_USAGE_score, SPLICE_JUNCTIONS_score, dnase_success, dnase_error, clndn, review_status, present, gnomad_variant_id, gnomad_rsids, gnomad_genome_{af,ac,an,ac_hom}, gnomad_exome_{af,ac,an,ac_hom}, error, pubmed_hits`

### Tier definitions

- **tier=1** (4 variants): the 4 wet-lab-testable Tier-1 candidates (rs801806, rs4293437, rs801809, rs2847163)
- **tier=2** (46 variants): top-50 ranked minus top-4
- **tier=3** (1,560 variants): all other SCN1A VUS in ClinVar

### Note on Tier-1 ranking

The 4 tier-1 candidates are at vus_rank 14, 17, 29, 30 — not at the very top of the splice ranking. Tier-1 selection was based on **mechanism + orthogonal evidence** (chromatin similarity to known pathogenic + brain RNA-seq effect + splice-site proximity), not on raw splice score alone. The top of the splice ranking is dominated by deep-intronic variants from `vus_tier2_candidates.csv` (see also `outputs/vus_tier2_candidates.csv` in this repo).