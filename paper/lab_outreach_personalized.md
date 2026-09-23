# SCN1A Lab Outreach — Personalized Emails (Ready to Send)

Generated 2026-09-22 with the latest DNase concentration finding (Exp 004).

Three emails, each ~250-350 words, ready to copy into your email. Contact info
and names need to be filled in by you before sending.

---

## EMAIL 1: Dr. Glen Saunders Carvill (Northwestern)

**To:** glenn.carvill at northwestern.edu (verify before sending)
**Subject:** SCN1A VUS re-scoring with AlphaGenome — DNase concentration signature + 4 Tier-1 candidates for minigene validation

Dr. Carvill,

I'm writing to share results from a systematic re-analysis of SCN1A variants using AlphaGenome (DeepMind's 2026 multimodal DNA model, https://alphagenome.google). Your 2018 *AJHG* paper on SCN1A non-coding variants and the 20N poison exon is directly relevant.

Two findings I'd like your lab's input on:

**(1) Pathogenic splicing-region variants in *SCN1A* show a chromatin-concentration signature beyond what AlphaGenome's splicing scores capture.** Our benchmark (n=60+87 variants across DMD/SCN1A/CFTR, Mann-Whitney U p=7.99 × 10⁻⁷, rank-biserial r=+0.467) shows pathogenic variants have *more concentrated* DNase-predicted chromatin perturbation at the variant position than benign intronic controls. This suggests AlphaGenome's chromatin predictors see something at pathogenic splice sites that generic splicing scores don't.

**(2) We identified 4 Tier-1 VUS candidates** (rs801806, rs4293437, rs801809, rs2847163) with explicit `splice_acceptor_variant` / `splice_donor_variant` annotations, ultra-rare in gnomAD (3 not observed, 1 at 0.003% AF), no prior publications, and high AlphaGenome splice scores (1.14, 1.06, 0.92, 0.92). These are direct candidates for minigene validation per the Sparber et al. 2023 protocol.

The full ranked candidate list (CSV), per-variant AlphaGenome scores, and analysis pipeline are at github.com/rollroyces/alphagenome-scn1a-pipeline (public; ~27 KB preprint draft also available).

Three questions:
1. Are any of these Tier-1 candidates known to your lab or your patient cohort?
2. Would your lab be open to applying the Sparber minigene assay to validate the highest-priority candidates?
3. Would you be open to co-authorship on a methods paper comparing AlphaGenome's chromatin concentration finding to your functional validation results?

No funding requested. The goal is to find out whether collaboration makes sense before either side invests more time. Happy to set up a 30-minute call or share additional detail in writing.

Best regards,
Royce Chi-Kit Chan
Independent researcher, Hong Kong

---

## EMAIL 2: Dr. P. Sparber / Corresponding author of Sparber 2023 (Research Center for Medical Genetics, Moscow)

**To:** corresponding author email (find via Springer paper page)
**Subject:** Extending your SCN1A minigene validation pipeline to 4 high-confidence AlphaGenome candidates

Dear Dr. Sparber,

Your 2023 *Human Genetics* paper validating 18 deep intronic *SCN1A* variants with the minigene assay is the most authoritative wet-lab pipeline for the exact question we've been working on computationally. I'm writing to propose extending your protocol to a new set of high-confidence candidates.

What we've done:
- Applied AlphaGenome (DeepMind's 2026 multimodal DNA model) to re-score all 1,610 SNV VUS in *SCN1A* in ClinVar
- Top 4 Tier-1 candidates (rs801806, rs4293437, rs801809, rs2847163) are splice-acceptor/donor variants, ultra-rare in gnomAD (3 not observed), no prior publications, with high AlphaGenome splice scores (1.14, 1.06, 0.92, 0.92)
- These candidates are in the same variant class as the 18 you validated in 2023

Why this matters:
- Your protocol is the gold-standard functional validation; computational predictions alone don't get into clinical practice
- Validating our AlphaGenome predictions against your experimental results would strengthen both our methods paper and any clinical interpretation tool built on AlphaGenome
- If our predictions are correct, ~5 of the top-10 VUS could be reclassified as likely pathogenic after experimental validation

What I'd like to discuss:
1. Would your lab be open to extending the 2023 minigene protocol to our 4 Tier-1 candidates?
2. Would you be open to co-authorship on a validation paper comparing AlphaGenome predictions to your experimental results?
3. Are there unpublished SCN1A VUS from your cohort that AlphaGenome could help prioritize?

The ranked candidate list and pipeline code are at github.com/rollroyces/alphagenome-scn1a-pipeline. Methods preprint (~5,000 words, ~27 KB) attached separately.

Happy to set up a 30-minute call or share additional detail in writing.

Best regards,
Royce Chi-Kit Chan
Independent researcher, Hong Kong

---

## EMAIL 3: Dr. Ingo Helbig (Children's Hospital of Philadelphia)

**To:** helbig at chop.edu (verify before sending)
**Subject:** SCN1A VUS re-scoring with AlphaGenome — integration with Epi25 unsolved Dravet cohort

Dr. Helbig,

Your 2021 *PLOS Genetics* review noted that ~20% of Dravet cases have no identified coding *SCN1A* variant. We've been systematically scoring *SCN1A* non-coding variants with AlphaGenome (DeepMind's 2026 multimodal DNA model) and would like to discuss integration with your Epi25 cohort work.

What we've found:
- 1,610 SNV VUS in *SCN1A* re-ranked by AlphaGenome splicing score
- 4 Tier-1 candidates (rs801806, rs4293437, rs801809, rs2847163) with explicit splice-acceptor/donor annotations, ultra-rare in gnomAD, no prior publications
- Cross-disease generalization: AUPRC > 0.98 for *SCN1A*, *SCN2A*, *MECP2*, *CFTR*, *DMD* on ClinVar pathogenic vs benign intronic SNVs
- New finding: pathogenic splice-region variants have more concentrated DNase chromatin perturbation at the variant position (n=60+87, Mann-Whitney U p=7.99 × 10⁻⁷, r=+0.467)

What I'd like to discuss:
1. Are any of these Tier-1 candidates already seen in Epi25 unsolved Dravet cases?
2. Could AlphaGenome splicing scores be used as PP3/BP4 computational evidence (ACMG/AMP 2015 standards) for reclassifying SCN1A VUS in your cohort?
3. Would you be open to co-authorship on a benchmark paper applying AlphaGenome to a broader set of developmental and epileptic encephalopathy (DEE) genes?

The ranked candidate list, full pipeline code, and methods preprint (~27 KB) are at github.com/rollroyces/alphagenome-scn1a-pipeline.

No funding requested. Happy to set up a 30-minute call or share additional detail in writing.

Best regards,
Royce Chi-Kit Chan
Independent researcher, Hong Kong

---

## Notes before sending

1. **Verify contact info.** The email addresses above are educated guesses from public sources. Confirm via the corresponding-author's recent paper before sending.
2. **Customize for your contact style.** These are formal/scientific. Adjust if you have a more direct relationship.
3. **One at a time.** Send to Carvill first (highest-fit lab), wait 1 week for reply, then Sparber, then Helbig.
4. **Include the preprint as a PDF attachment** (or link to GitHub repo). The preprint is at `paper/preprint.md` — convert to PDF first.
6. **Follow up.** If no reply in 2 weeks, send a single polite follow-up email. If still no reply, move on.