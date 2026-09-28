---
title: SCN1A VUS Re-Scoring
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: streamlit
sdk_version: 1.64.0
app_file: app.py
pinned: false
license: mit
short_description: Browse 1,610 SCN1A VUS scored with AlphaGenome + DNASE
---

# SCN1A VUS Re-Scoring — Clinician Web UI

<!-- DEPLOYED_EMBED_START -->
<iframe
  src="https://huggingface.co/spaces/RROL/alphagenome-scn1a-app"
  width="100%" height="900" frameborder="0" allow="clipboard-read; clipboard-write"
></iframe>
<!-- DEPLOYED_EMBED_END -->

A clinician-facing web UI for the **1,610 SCN1A variants of uncertain significance (VUS)**
re-scored with [AlphaGenome](https://alphagenome.google) (DeepMind, 2025) splicing + DNASE scorers
and cross-referenced against gnomAD v4.1 and PubMed.

The dataset is published on Hugging Face:
**[RROL/scn1a-vus-alphagenome](https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome)**
(CC-BY-4.0, 30 columns).

## Pages

| # | Page | What it does |
|---|------|--------------|
| 1 | **Browse Dataset** | Filter all 1,610 VUS by consequence, AlphaGenome scores, gnomAD AF, and tier. Click a row for full variant detail. |
| 2 | **Tier-1 Candidates** | Detailed evidence cards for the **4 wet-lab-ready splice-region VUS** prioritized for the Carvill lab minigene validation (rs801806, rs4293437, rs2847163, rs801809). |
| 3 | **Tier-2 Candidates** | Detailed evidence cards for the **13 curated Tier-2 candidates** spanning 3 mechanism classes (deep intronic cryptic splice, dual-mechanism missense, non-coding transcript overlap). All absent from gnomAD, 0 SCN1A-specific PubMed hits. |
| 4 | **Submit VCF** | Paste or upload a VCF; score each variant with the live AlphaGenome API. **Uses real API quota** — the API key you enter is held in session state only, never persisted. |

The app works without an API key for Pages 1–3. Page 4 requires a free key from
[https://alphagenome.google/api](https://alphagenome.google/api).

## Methodology

Full pipeline, raw score data, per-experiment README, and lab outreach templates are in the
GitHub repository:

**[github.com/rollroyces/alphagenome-scn1a-pipeline](https://github.com/rollroyces/alphagenome-scn1a-pipeline)**

Key references inside the repo:

- `outputs/tier1_actionable_features.md` — wet-lab brief for the 4 Tier-1 candidates
- `outputs/vus_tier2_candidates.md` — 13 Tier-2 mechanism-stratified candidates
- `paper/preprint.md` — full methods draft (~5,000 words)
- `research_notebook/experiments/013_tier1_actionable_features/` — Exp 013 methodology
- `research_notebook/experiments/010_tier2_candidates/` — Exp 010 methodology

## Tier definitions

- **tier = 1** (4 variants): the 4 wet-lab-testable Tier-1 candidates (rs801806, rs4293437, rs801809, rs2847163)
- **tier = 2** (46 variants): top-50 splice-ranked minus top-4
- **tier = 3** (1,560 variants): all other SCN1A VUS in ClinVar

Tier-1 selection is based on **mechanism + orthogonal evidence** (chromatin similarity to known
pathogenic + brain RNA-seq effect + splice-site proximity), **not** on raw splice score alone.
The 13 Tier-2 cards on Page 3 are a curated subset of the 46 tier-2 variants, stratified into
3 mechanism classes (5 deep intronic, 5 dual-mechanism missense, 3 non-coding transcript).

## Running locally

```bash
git clone https://github.com/rollroyces/alphagenome-scn1a-pipeline
cd alphagenome-scn1a-pipeline
pip install -r requirements.txt
streamlit run app.py
```

The app launches on `http://localhost:8501`. Page 4 (Submit VCF) additionally needs
`pip install alphagenome` and a free API key.

## Free public hosting (no HF Pro required)

Hugging Face Spaces for Streamlit require a **Pro subscription** as of 2026
(see [HF Spaces free-tier changes](https://toolfreebie.com/hugging-face-spaces-free-gpu/)).
The free alternatives that work with this codebase as-is:

- **[Streamlit Community Cloud](https://share.streamlit.io)** — GitHub-connected free hosting;
  connect this repo, point at `app.py`, deploy in ~2 minutes. **This is the recommended
  hosting target.** Unlimited public apps, free tier, sleeps on inactivity.
- **Render / Railway / Fly.io** — generic Docker hosts; ship the repo as a container with
  `streamlit run app.py --server.port $PORT --server.address 0.0.0.0` as the entrypoint.

This Space repo (`RROL/alphagenome-scn1a-app`) stores the deployable artifacts but cannot
run as a live Streamlit Space on the RROL free-tier account. The Space runtime is paused
(`errorMessage: Quota exceeded for flavor cpu-basic`).

## Citation

```
Royce Lam. (2026). SCN1A VUS Re-Scoring with AlphaGenome (Revision 8db4888).
Hugging Face Spaces. https://huggingface.co/spaces/RROL/alphagenome-scn1a-app
```

If you use this work, please cite the underlying dataset:

```
Royce Lam. (2026). SCN1A VUS Re-Scoring with AlphaGenome (Revision 8db4888). Hugging Face.
https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome
```

## License

- App code: MIT
- Dataset: CC-BY-4.0
- AlphaGenome SDK & models: see [DeepMind AlphaGenome terms](https://alphagenome.google)