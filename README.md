# AlphaGenome work — rare disease non-coding variant scoring

Re-scoring ClinVar variants of uncertain significance (VUS) in *SCN1A* with [AlphaGenome](https://alphagenome.google) (DeepMind, 2025) to identify novel candidate pathogenic variants for the unsolved fraction of Dravet syndrome.

**Status:** v0.1 — methods paper draft, candidate list ready for outreach.

## Headline results

| Metric | Value |
|--------|-------|
| AlphaGenome benchmark AUPRC on SCN1A | **0.9833** |
| Top-5% precision | **100%** |
| AUROC | 0.9896 |
| Variants scored | 591 (216 pathogenic + 375 benign controls) |
| API call success rate | 591/591 (100%) |
| VUS re-scored | 1,610 |
| **Novel high-impact candidates** | **61** (3.8% of VUS) |
| Top-20 candidates with prior literature | **0 / 20** |

**The core finding:** 61 SCN1A VUS are flagged as high-impact by AlphaGenome, none have prior SCN1A publications, the 8 that appear in gnomAD are absent or ultra-rare, and 19 are missense variants predicted to cause dual-mechanism (coding + splice) pathogenicity.

See `paper/preprint.md` for the full draft (~5,000 words), `paper/candidate_report.md` for the lab-ready variant report, and `paper/outreach_template.md` for email templates to 3 target labs (Carvill, Sparber, Helbig).

## Quickstart

1. Get a free AlphaGenome API key at https://alphagenome.google/api.
2. `echo 'YOUR_KEY' > .alphagenome_key && chmod 600 .alphagenome_key` (in project root).
3. `uv venv --python 3.13 .venv && uv pip install --python .venv/bin/python alphagenome pandas numpy jupyter matplotlib seaborn scikit-learn requests tqdm`
4. `bash scripts/_run_with_key.sh scripts/benchmark_scn1a_live_api.py` (benchmark, ~7 min)
5. `bash scripts/_run_with_key.sh scripts/rescore_vus.py` (VUS re-scoring, ~18 min)

The full project documentation below was written for the first-pass walk-through; the GitHub-facing summary above is the canonical entry point.

---

## Honest framing

This is **not** a Nobel Prize project. It is a rigorous computational study
that can:
1. Produce a credible preprint within 90 days
2. Become a peer-reviewed paper within 6–9 months
3. Serve as a portfolio / proof-of-work for finding a biology collaborator
4. Become a foundation for an actual research program (10+ year horizon)

The Nobel Prize in Physiology/Medicine rewards discoveries, not infrastructure.
This work could become infrastructure *for* a discovery — that is the realistic
path.

## Project layout

```
alphagenome-work/
├── .venv/                    # Python 3.13 venv (uv-managed)
├── tutorials/                # Official DeepMind tutorial notebooks
│   └── alphagenome-upstream/ # git clone of google-deepmind/alphagenome
├── data/                     # Reference data (gnomAD, ClinVar, GENCODE)
├── notes/                    # Learning log, paper notes, biology primer
├── figures/                  # Plots, locus diagrams, calibration curves
├── outputs/                  # Pipeline outputs (ranked candidate lists)
└── README.md                 # This file
```

## Phasing

### Phase 1 (Days 1–10): Tutorials + biology primer
- Run all 10 official tutorial notebooks
- Document what AlphaGenome can / can't do
- Write a primer on GWAS / VCF / Mendelian disease for the future-us

### Phase 2 (Days 11–20): Disease selection
- Survey 5–10 candidate unsolved rare diseases
- Pick one with: public WGS data, established unsolved fraction, clean labels

### Phase 3 (Days 21–60): Pipeline build
- VCF → AlphaGenome API + Atlas → ClinVar annotation → ranked candidates
- Calibration analysis (isotonic / Platt)
- Precision/recall per modality

### Phase 4 (Days 61–80): Real analysis + paper draft
- Apply to unsolved fraction of chosen disease
- Generate candidate list
- Draft preprint

### Phase 5 (Days 81–90): Outreach
- Identify labs working on the disease
- Send ranked candidate list
- Find collaborator

## Tools

- Python 3.13 (uv-managed venv)
- `alphagenome` SDK (gRPC API client)
- pandas / numpy / scikit-learn for analysis
- jupyter for exploration
- matplotlib / seaborn for figures
- ClinVar / gnomAD / GENCODE for ground truth

## What I do vs what Hermes does

**You do:**
- Run the tutorials
- Read the papers
- Write the pipeline code
- Run the analysis
- Draft the paper

**Hermes helps with:**
- Pipeline design review
- Biology translation (papers → concepts)
- Stats debugging
- Paper structure / prose polish
- Outreach email drafts

Hermes does NOT do the work for you. You will write code, run it, read results.

## First concrete action

1. Get an API key at https://alphagenome.google/api (free, non-commercial)
2. Run `quick_start.ipynb` end to end
3. Report back: what worked, what confused you, what surprised you

That tells us how to phase Phase 2.
