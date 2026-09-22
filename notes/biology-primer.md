# Biology primer for the data scientist

A cheat sheet I wish I'd had when starting. Written for someone who knows
data science cold but has zero biology background.

## The central dogma (in 5 lines)

```
DNA  →  RNA  →  Protein
```

- **DNA**: the long-term storage. 4 letters: A, T, G, C. Two complementary strands wound into a double helix. Human genome is ~3.2 billion letters.
- **RNA**: a working copy of a gene. Also 4 letters (A, U, G, C — U replaces T). Single-stranded.
- **Protein**: the machine that does work. 20 amino acids.

When a gene is "expressed," the cell copies DNA→RNA (transcription), then
RNA→protein (translation). "Gene regulation" controls *when*, *where*, and
*how much* this happens.

## The genome has two parts

1. **Coding regions** (~2% of the genome): directly specify protein sequence. Made of **exons** separated by **introns**.
2. **Non-coding regions** (~98%): everything else. Used to be called "junk DNA" — now known to contain:
   - **Promoters**: where transcription starts
   - **Enhancers**: amplify expression of nearby (or far) genes
   - **Silencers**: suppress expression
   - **Splice sites**: tell the cell where to cut introns out of RNA
   - **Transcription factor binding sites**: where regulatory proteins attach
   - **Non-coding RNAs**: functional RNA molecules that aren't translated

**This is what AlphaGenome is good at predicting: what happens in the 98%.**

## A single nucleotide variant (SNV) — what does it actually do?

When a single base pair changes (e.g., A→G at position chr7:117,559,590), the consequence depends on *where* it lands:

| Location | Possible effect |
|---|---|
| Coding exon (synonymous) | Same amino acid, usually benign |
| Coding exon (missense) | Different amino acid, may affect protein |
| Coding exon (nonsense) | Premature stop codon, truncated protein |
| Splice site (donor/acceptor ±1,2) | Exon skipping, intron retention — usually pathogenic |
| Enhancer / promoter | Alters expression level — clinically important, hard to detect |
| TF binding site | Disrupts regulatory protein binding — same as above |
| Deep intronic | Often no effect, occasionally creates cryptic splice site |

**AlphaGenome predicts the regulatory consequences — the bottom three rows.**
**ClinVar catalogs known pathogenic variants with clinical evidence.**

## Variant nomenclature

`chr7:117,559,590 G>A` means:
- chromosome 7
- position 117,559,590
- reference allele G
- observed allele A

`NM_001234.5:c.76A>C` means:
- transcript NM_001234.5 (GenBank ID, version 5)
- at coding position 76
- A→C

You'll see both styles. Stick with chr:pos:ref>alt (the VCF-style) for
AlphaGenome work.

## VCF — Variant Call Format

Standard format for storing variants. Columns:

```
#CHROM  POS     ID   REF  ALT  QUAL  FILTER  INFO
chr7    117559590  .  G    A    99    PASS    .
```

Each row = one variant. For AlphaGenome, you only need CHROM, POS, REF, ALT.

## How rare disease genetics actually works

A "rare disease" = a Mendelian disorder, usually caused by a single gene
breaking. Examples:
- Cystic fibrosis (CFTR)
- Huntington's (HTT)
- Duchenne muscular dystrophy (DMD)
- ~7,000 known rare diseases total

When a patient presents with symptoms suggesting a genetic cause, the
clinical workflow is:
1. Whole exome sequencing (WES) — reads coding regions only
2. If unsolved: whole genome sequencing (WGS) — reads everything
3. Filter variants: rare (gnomAD AF < 0.01), in known disease gene, predicted damaging
4. Classify per ACMG criteria: pathogenic / likely pathogenic / VUS / likely benign / benign

**~50–70% of WGS-negative rare disease patients remain unsolved.**
Hypothesis: many of these unsolved cases have causal variants in
non-coding regions that current tools can't interpret.

**This is where AlphaGenome fits.**

## What AlphaGenome actually predicts

Given a 1 Mb DNA sequence, it outputs numerical tracks for 11 modalities
across many cell types / tissues:

- **RNA-seq** (gene expression level)
- **CAGE / PRO-cap** (transcription start sites)
- **Splice sites, splice usage, splice junctions** (how RNA is processed)
- **DNase, ATAC-seq** (chromatin accessibility)
- **Histone marks** (H3K4me3, H3K27ac, etc.)
- **Transcription factor binding** (where TFs land)
- **Chromatin contact maps** (3D structure)

For variant scoring: you submit a variant, the model computes ref and alt
predictions, and the *delta* tells you how much that single letter change
disrupts the predicted biology.

## What AlphaGenome does NOT do

- Predict protein structure (use AlphaFold)
- Predict mRNA stability / translation efficiency / codon optimality
- Predict cell-to-cell variation (it gives bulk tissue signals)
- Predict effects of structural variants (only SNVs and small indels)
- Tell you if a variant is "pathogenic" — it tells you if a variant is
  *functionally disruptive*; clinical interpretation is a separate step

## Why the field is hard

For most variants, you can't just say "AlphaGenome says it's bad → it's
pathogenic." You need:
- Population frequency (gnomAD)
- Gene-level constraint (pLI, LOEUF)
- Clinical case reports (ClinVar, HGMD, case series)
- Functional assays (in vitro / in vivo validation)

A high AlphaGenome score is one input among many. The model is most useful
for *prioritizing* variants worth follow-up, not for definitive calls.

## Key public databases (you'll use these)

| DB | What's in it | Where |
|---|---|---|
| gnomAD | Population allele frequencies, ~800k exomes+genomes | https://gnomad.broadinstitute.org/ |
| ClinVar | Clinical variant classifications (pathogenic / benign) | https://www.ncbi.nlm.nih.gov/clinvar/ |
| GENCODE | Gene annotation (where exons/introns/enhancers are) | https://www.gencodegenes.org/ |
| GTEx | Tissue-specific gene expression | https://gtexportal.org/ |
| ENCODE | Functional genomics tracks (histone marks, TF binding) | https://www.encodeproject.org/ |
| UCSC Genome Browser | Browser for visualizing all the above | https://genome.ucsc.edu/ |
| SpliceAI | Specialized splicing effect predictor | Illumina BaseSpace / standalone |
| CADD / REVEL / AlphaMissense | Other variant effect predictors | various |

## Reading list (start here, in this order)

1. **Avsec et al. 2026** (Nature) — AlphaGenome paper
2. **100,000 Genomes Project main results** — what WGS finds in rare disease
3. **Claussnitzer et al. 2020** — review on non-coding variant interpretation
4. **ACMG/AMP 2015 standards** — how variants are clinically classified
5. **Cooper 2024 review** — current state of rare disease genomics

I'll point you to specific papers as we go.

## Quick glossary

- **Mendelian**: caused by a single gene. Dominant = one copy enough; recessive = need two.
- **Penetrance**: fraction of carriers who actually develop disease.
- **Expressivity**: how severe the disease is when it manifests.
- **Haploinsufficiency**: one working copy not enough, causes disease.
- **VUS**: variant of uncertain significance — ClinVar's "we don't know" bucket.
- **Compound het**: two different pathogenic variants in the same gene, one on each copy (recessive pattern).
- **De novo**: variant that arose in the patient, not inherited from parents.
