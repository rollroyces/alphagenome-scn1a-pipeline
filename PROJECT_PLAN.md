# SCN1A / Dravet syndrome — project plan

## Disease background

- **Dravet syndrome**: severe developmental and epileptic encephalopathy (DEE), onset in infancy, refractory seizures, developmental regression. Incidence ~1:15,000.[1]
- **Genetic basis**: ~80% of clinically diagnosed Dravet patients carry pathogenic loss-of-function variants in **SCN1A** (voltage-gated sodium channel, Nav1.1, GABAergic interneurons).[1]
- **Unsolved fraction**: ~20% of clinically diagnosed Dravet patients have no SCN1A coding variant identified even after WES/WGS.[2]

## Why this is the right disease for our project

1. **Single gene** — exhaustive scope, low scope-creep risk
2. **Mechanism established** — noncoding variants activate poison exons → aberrant splicing → haploinsufficiency (Carvill 2018, Voskobiynyk 2020 mouse model, Sparber 2023)
3. **AlphaGenome's strongest modality** — splicing prediction, the exact variant class we're hunting
4. **Validation path exists** — Sparber 2023 minigene assay validated for all SCN1A poison exons
5. **Active research community** — multiple labs with WGS cohorts and wet-lab capacity

## Specific scientific anchor — the 20N poison exon

The **20N poison exon** (also called exon 20N) is the best-characterized case:
- Located in intron 20 of SCN1A
- Flanked by an SRSF1 binding motif that normally suppresses inclusion
- Variants disrupting SRSF1 binding → aberrant inclusion of 20N → premature stop codon → NMD → SCN1A haploinsufficiency
- Carvill 2018: **7 noncoding variants** identified in Dravet patients, **5 of which** activate 20N inclusion
- Sparber 2023: validated minigene assay covering 20N and other poison exons, tested 18 deep-intronic SNVs

**This is our positive-control set.** If AlphaGenome correctly predicts high splice-disruption scores for the Carvill 2018 variants, the model is working as expected on real pathogenic variants.

## Project phases

### Phase 2 (Days 11–30): Data acquisition + pipeline scaffold

**Week 1: SCN1A locus + ClinVar**
- [ ] Download SCN1A gene interval (hg38 chr2:165,984,640–166,149,114)
- [ ] Pull all SCN1A variants from ClinVar (full release + REST API subset)
- [ ] Pull all SCN1A variants from gnomAD (population frequencies)
- [ ] Catalog all SCN1A exons / introns / known poison exons from GENCODE + literature
- [ ] Catalog all known pathogenic SCN1A variants from Carvill 2018, Sparber 2023, HGMD

**Week 2: Atlas pipeline**
- [ ] Use `AtlasClient.get_interval()` to score ALL ~150,000 possible SNVs in SCN1A locus with all splicing scorers
- [ ] Cache results to `data/atlas_scn1a.h5ad`
- [ ] Compute summary statistics (score distributions by region: exon / intron / poison-exon-flanking)

**Week 3: Live API validation**
- [ ] Take the 7 Carvill 2018 variants and the c.4853-25 T>A variant — submit each to `score_variant` with all 18 recommended scorers
- [ ] Compare live predictions against Atlas pre-computed scores
- [ ] Verify AlphaGenome correctly ranks the known pathogenic variants as high-impact

**Deliverable:** `outputs/scn1a_atlas_scores.tsv` — every SNV in SCN1A scored for splicing impact, ready for candidate ranking.

### Phase 3 (Days 31–50): Benchmark + validation

**Week 4: Benchmark against known pathogenic variants**
- [ ] Build a labeled set:
  - Positive: 7 Carvill 2018 variants + 18 Sparber 2023 variants + c.4853-25 T>A = ~26 known pathogenic splicing variants
  - Negative: common (gnomAD AF > 0.01) intronic variants in SCN1A = ~thousands of benign controls
- [ ] Compute AUPRC per modality (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS)
- [ ] Compute calibration curves (do raw scores correspond to probability of pathogenicity?)
- [ ] Compare AlphaGenome's splicing performance against SpliceAI on the same labeled set

**Week 5: Novel candidate identification**
- [ ] Rank all SCN1A intronic variants by AlphaGenome splicing impact
- [ ] Filter: rare in gnomAD (AF < 0.001), high splicing impact, near known poison exons or splice sites
- [ ] Cross-reference with literature for any "VUS" or "novel" annotations
- [ ] Produce top-20 candidate list with predicted mechanism

**Deliverable:**
- `outputs/benchmark_results.json` — precision/recall/AUPRC per modality
- `outputs/scn1a_top_candidates.csv` — top 20 candidates with mechanism prediction
- `figures/scn1a_locus_plot.png` — locus diagram showing score distribution across SCN1A
- `figures/calibration_curves.png` — calibration plots

### Phase 4 (Days 51–75): Paper draft + figures

- [ ] Draft preprint on bioRxiv format
- [ ] Figures: locus plot, calibration, benchmark table, top-20 candidate table
- [ ] Methods: pipeline reproducibility, version pinning, data sources
- [ ] Discussion: comparison to SpliceAI, limitations, future directions

**Deliverable:** `paper/preprint.md` — full draft ready for bioRxiv submission.

### Phase 5 (Days 76–90): Outreach

- [ ] Identify 3–5 labs working on SCN1A / Dravet / noncoding epilepsy variants
- [ ] Draft personalized outreach emails with the candidate list
- [ ] Submit to bioRxiv
- [ ] Follow up with interested labs to establish collaboration for validation

## Specific candidate labs to contact (initial list)

| Lab | PI | Why |
|---|---|---|
| Helbig lab (CHOP) | Ingo Helbig | PLOS Genetics 2021 review on poison exons, active in SCN1A noncoding variant interpretation |
| Carvill lab (Lurie Children's) | Gemma Carvill | Lead author on Carvill 2018 (original poison exon paper), maintains SCN1A noncoding variant database |
| Sparber lab (Research Centre for Medical Genetics) | Peter Sparber | Lead author on Sparber 2023 (minigene assay for SCN1A poison exons) |
| Marini lab (Meyer Children's Hospital) | Carla Marini | Italian Dravet cohort, multiple SCN1A variant publications |
| Epi25 collaborative | multi-PI | Large-scale WGS in epilepsy, has unsolved Dravet cases |

We contact them in order of likelihood to engage.

## What this can realistically produce (revised estimate)

| Outcome | Probability | Deliverable |
|---|---|---|
| Methods paper on bioRxiv | ~80% | "AlphaGenome splicing predictions benchmarked on SCN1A noncoding variants" |
| Peer-reviewed methods paper | ~50% | Same content, submitted to Genome Biology / Bioinformatics / NAR |
| Novel candidate variant flagged | ~40% | One of our top-20 hits turns out to be a true novel SCN1A noncoding pathogenic variant (no wet-lab validation, just prediction) |
| Wet-lab collaboration established | ~30% | A lab agrees to validate one of our top candidates via minigene assay |
| Functional validation of a candidate | ~15% | A candidate variant is experimentally confirmed |
| Reclassification of an existing VUS | ~10% | One of our high-scoring variants gets reclassified from VUS to likely pathogenic |

The probability of any *single* breakthrough is modest. The probability of *some* useful outcome is high. That's the honest math.

## Sources

[1] Helbig I, Goldberg E (2021) "The dose makes the poison—Novel insights into Dravet syndrome and SCN1A regulation through nonproductive splicing" PLOS Genet 17(1):e1009214. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009214

[2] Carvill GL et al. (2018) "Aberrant Inclusion of a Poison Exon Causes Dravet Syndrome and Related SCN1A-Associated Genetic Epilepsies" AJHG 103(6):1022-1029. https://www.sciencedirect.com/science/article/pii/S0002929718303999

[3] Sparber P et al. (2023) "Functional investigation of SCN1A deep-intronic variants activating poison exons inclusion" Hum Genet 142:1043-1053. https://link.springer.com/article/10.1007/s00439-023-02564-y

[4] Zhang G et al. (2025) "Dravet syndrome: novel insights into SCN1A-mediated epileptic neurodevelopmental disorders within the molecular diagnostic-therapeutic framework" Front Neurosci 19:1634718. https://pmc.ncbi.nlm.nih.gov/articles/PMC12326748/
