# paper/preprint.md — Phase 2 polish checklist

Items to address before bioRxiv submission, in priority order:

## A. Author & affiliation (must do)
- [ ] Real affiliation for Royce
- [ ] Real ORCID
- [ ] Either remove Hermes Agent author OR clarify AI-tool usage per target journal's policy
- [ ] Email for correspondence

## B. Methods clarifications (must do)
- [ ] Section 2.4: be explicit that direct SpliceAI head-to-head was infeasible; cite AlphaGenome's published SpliceAI comparison more directly
- [ ] Section 2.5: clarify "RCS writes pipeline + drafted paper" vs "Hermes drafts initial draft, Royce verifies and edits"
- [ ] Add Gencode version, ClinVar release date, and download date for full reproducibility

## C. Results additions (nice to have)
- [ ] Add Section 3.6 cross-disease benchmark table (already there — verify formatting)
- [ ] Add Section 3.7 referencing research_notebook/experiments/001 ISM negative result as a future direction
- [ ] Add Limitations: training-set leakage explicitly not controlled

## D. Figures
- [ ] Verify all figures render properly
- [ ] Add cross-disease AUPRC bar chart reference
- [ ] Add AUPRC benchmark figures reference

## E. References
- [ ] Verify all DOIs are correct
- [ ] Add AlphaGenome paper citation (Avsec et al. 2026, Nature) with full author list
- [ ] Add ClinVar / gnomAD URLs

## F. Submission prep
- [ ] Decide on target journal: bioRxiv (preprint) → Genome Biology / Bioinformatics / Nature Genetics?
- [ ] bioRxiv needs license declaration (CC-BY-4.0 or CC-BY-NC-4.0)
- [ ] Need subject area tags

## G. Defer / not blocking
- SpliceAI comparison (would require cloud GPU)
- Wet-lab validation
- Other modalities (in research_notebook/experiments/002_ if that runs)
