# SCN1A VUS re-scoring — outreach email draft

Use this template once the candidate list is finalized. Each target lab gets a customized version.

---

## Target labs (priority order)

### 1. Carvill lab — Glen Saunders Carvill (lead author, Carvill 2018)
- **Institution:** Northwestern University Feinberg School of Medicine
- **Why:** Identified the original SCN1A poison exon 20N; their lab has the functional pipeline and patient cohort
- **Hook:** "We've applied AlphaGenome to the 1,610 SNV VUS in SCN1A and identified a short list of high-confidence candidates for functional validation. Top candidates are deep intronic variants predicted to disrupt the same poison exon pathway your lab characterized in 2018."
- **What we offer:** ranked candidate list (CSV), pipeline code, AlphaGenome scores for every VUS
- **What we ask:** access to your SCN1A functional validation pipeline (minigene assay per Sparber 2023), co-authorship on validation paper
- **Find contact:** https://www.feinberg.northwestern.edu/faculty-profiles/az/profile-20989.html

### 2. Sparber lab — P. Sparber (lead author, Sparber 2023)
- **Institution:** Research Center for Medical Genetics, Moscow
- **Why:** Validated 18 deep intronic SCN1A variants experimentally; the most authoritative wet-lab pipeline for this exact question
- **Hook:** "We applied AlphaGenome to re-score all 1,610 SNV VUS in SCN1A. The top 10 candidates overlap with the deep-intronic variant classes you validated in 2023. We'd like to extend your protocol to the remaining high-priority VUS."
- **What we offer:** ranked candidate list, validation-ready VCF subset
- **What we ask:** access to minigene assay, protocol exchange, co-authorship
- **Find contact:** https://link.springer.com/article/10.1007/s00439-023-02564-y (corresponding author in paper)

### 3. Helbig lab — Ingo Helbig (senior author, Helbig 2021 review)
- **Institution:** Children's Hospital of Philadelphia (CHOP)
- **Why:** Senior voice in Dravet syndrome genetics; sees thousands of unsolved patient cases through Epi25
- **Hook:** "We're systematically scoring SCN1A VUS with AlphaGenome. We have ~50 high-confidence candidates for functional validation that could immediately be checked against your Epi25 cohort."
- **What we offer:** ranked candidate list, integration plan with Epi25
- **What we ask:** clinical context for high-impact VUS, co-authorship on benchmark paper
- **Find contact:** https://www.chop.edu/doctors/helbig-ingo

---

## Email template (customize per lab)

```
Subject: SCN1A VUS re-scoring with AlphaGenome — high-confidence splicing candidates for validation

Dear Dr. [Last Name],

I'm writing to share results from a systematic re-analysis of SCN1A variants
using AlphaGenome (DeepMind's multimodal DNA model for non-coding variant
interpretation, https://alphagenome.google). Your [lab's 2018 SCN1A poison exon
work / 2023 minigene validation study / 2021 review on Dravet unsolved cases]
is directly relevant to what we've found.

Summary of work:
- AlphaGenome's SPLICE_SITES variant scorer discriminates pathogenic from
  benign intronic SCN1A variants with AUPRC = 0.9833 and top-5% precision =
  100% on a 591-variant ClinVar benchmark.
- We applied the same pipeline to all 1,610 SNV VUS currently in ClinVar for
  SCN1A. The top [N] candidates are deep intronic variants predicted to
  disrupt splicing — directly relevant to [the poison exon pathway you
  characterized / the variants you validated / the unsolved Dravet cases you
  reviewed].
- The full ranked candidate list (CSV), per-variant AlphaGenome scores, and
  the analysis pipeline are attached / available at github.com/[...].

What I'd like to discuss:
1. Are [N] of these top candidates known to you / your lab / your patient
   cohort?
2. Would your lab be open to applying [your minigene assay / your patient
   cohort] to validate the highest-priority candidates?
3. Would you be open to co-authorship on a methods paper comparing AlphaGenome
   splicing predictions to [your experimental validation results / ClinVar
   gold-standard labels]?

No funding or commitment is requested at this stage. The candidate list and
pipeline are ready to share; the goal of this email is to find out whether
collaboration makes sense before either side invests more time.

Happy to set up a 30-minute call or to share additional methods detail in
writing.

Best regards,
Royce Chi-Kit Chan
Independent researcher, Hong Kong
+852 [phone] / [email]
```

---

## Customization checklist before sending

For each lab:
- [ ] Replace bracketed placeholders with lab-specific facts
- [ ] Verify the actual corresponding author and current institution
- [ ] Reference their specific most recent paper (not just Carvill 2018)
- [ ] Have a concrete next step (call, follow-up email)
- [ ] Sign with real contact info

## What NOT to do

- Do not promise AlphaGenome is "ground truth" — it's a prediction; the lab's
  wet-lab validation is the truth
- Do not offer anything you can't deliver (e.g. specific clinical interpretation
  advice; that's the lab's job)
- Do not send to more than 5 labs simultaneously — it dilutes any individual
  collaboration and looks mass-spam
