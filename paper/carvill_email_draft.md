# Carvill Lab Outreach — Final Draft Email

**To:** glenn.carvill@northwestern.edu (verify before sending — search "Glenn Carvill Northwestern Neurology" first)  
**CC:** www.carvill-lab.org if exists  
**Attachment:** `outputs/tier1_actionable_features.md` (the lab handout)

---

### Subject

SCN1A VUS — 2 ranked candidates ready for minigene validation (1-page attached)

---

### Email body

Dr. Carvill,

I've re-scored all 1,610 SCN1A VUS in ClinVar with AlphaGenome and ranked the top 17 by 3 orthogonal evidence types (chromatin similarity to known pathogenic patterns, brain RNA-seq effect, splice-site proximity). My top 2 candidates for your published Sparber 2023 minigene assay are **rs801806** (cosine 0.24 to known path patterns, splice acceptor disruption, exon 16) and **rs4293437** (9× stronger brain RNA-seq LFC, splice acceptor, exon 4).

The attached 1-page brief has the full ranking, mechanism hypotheses, AlphaGenome scores, and gnomAD data. None of the 4 Tier-1 candidates have published functional data.

Would a 15-minute Zoom work to discuss whether these are testable in your lab? I can also send the pre-print when it's posted to bioRxiv.

No funding requested. Co-authorship on any resulting publication would be appropriate.

Best regards,  
Royce Chan  
Independent researcher, Hong Kong  
github.com/rollroyces/alphagenome-scn1a-pipeline

---

# Why this version

| Section | Length | Purpose |
|---|---|---|
| Opening line | 1 sentence | Shows I've done work; names the paper they published |
| Method summary | 1 sentence | "Re-scored 1,610 VUS with AlphaGenome, ranked by 3 evidence types" |
| Specific ask | 1 sentence | "Top 2 candidates are rs801806 and rs4293437" |
| Evidence for the rank | — | Cosine, LFC, mechanism, exon location |
| Attachment | — | "Full 1-page brief attached" |
| Next step | 1 sentence | "15-min Zoom" |
| Funding/co-authorship | 1 sentence | "No funding requested, co-authorship appropriate" |

**Total: ~140 words in body. Reads in 30 seconds.**

# What I deliberately cut

- ✗ Long background on AlphaGenome (he already knows what it is)
- ✗ Multiple paragraphs about our 8-gene benchmark (irrelevant to him)
- ✗ The full 5-page candidate report (just attach the 1-page)
- ✗ The other 13 Tier-2 candidates (we can mention on the Zoom)
- ✗ The negative findings (Exp 003, 005, 006, 007 nulls — not interesting yet)
- ✗ Promises of a paper (paper doesn't help him)
- ✗ Gushing about his prior work (he knows his work is good)

# What he gets in the email

1. **3 sentences** he can read in 30 seconds
2. **2 specific variants** with concrete evidence
3. **A 1-page attachment** he can look at if interested
4. **A clear ask** (15-min Zoom)
5. **No obligation** (no funding, no co-author requirement)

# What to do before sending

1. **Verify his email** — search Google: "Glenn Carvill Northwestern email"
2. **Check if lab website exists** — for the CC line
3. **Verify the attachment opens** — `outputs/tier1_actionable_features.md` should be a clean readable markdown
5. **Get a bioRxiv link if possible** — bioRxiv credibility boost (skip if no time)

# Companion email for Sparber (optional, send if Carvill doesn't reply in 2 weeks)

```
Subject: SCN1A minigene candidates — extending your 2023 validation pipeline

Dr. Sparber,

Your 2023 minigene validation of 18 SCN1A deep intronic SNVs is the 
gold-standard pipeline for exactly the question we've been working on 
computationally.

We've identified 4 Tier-1 SCN1A VUS candidates with AlphaGenome that 
are in the same variant class as the 18 you validated:
- rs801806 (splice_acceptor_variant, exon 16)
- rs4293437 (splice_acceptor_variant, exon 4)  
- rs801809 (splice_donor_variant, exon 14)
- rs2847163 (splice_donor_variant, exon 14)

Would your lab be open to extending the 2023 protocol to these? 
1-page brief attached.

Best regards,
Royce Chan
Independent researcher, Hong Kong
```

# Companion email for Helbig (also optional)

Helbig is more clinical-genetics focused. He might not do minigene assays but might forward to a lab that does. Worth a separate short email.

```
Subject: SCN1A VUS re-classification candidates — 4 ranked for clinical review

Dr. Helbig,

I've identified 4 SCN1A VUS candidates (rs801806, rs4293437, rs801809, 
rs2847163) using AlphaGenome that have explicit splice_region_variant 
annotations and ultra-rare allele frequencies. These are candidates 
for re-classification pending functional validation.

Would your clinical genomics team have a workflow for triaging these 
for the Dravet variant curation effort? 1-page brief attached.

Best regards,
Royce Chan
```

# Final action plan

1. **Verify Carvill's email** (15 min)
2. **Send Carvill email** with attachment (5 min)
4. **If no reply in 2 weeks** → send Sparber email
5. **If no reply in 6 weeks** → send Helbig email

Total time: 20-25 minutes for the highest-impact actions in this entire project.