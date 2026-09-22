# Learning log

Day-by-day notes on what I learned, what surprised me, what to revisit.

## Day 1 (today)

### What I did
- Set up project at `~/projects/alphagenome-work/`
- Created Python 3.13 venv via `uv venv` (note: PYTHONPATH hijack quirk — must use `uv pip install --python .venv/bin/python`)
- Installed `alphagenome` SDK + pandas / numpy / scikit-learn / jupyter / matplotlib / seaborn
- Cloned official tutorial repo (10 colabs in `tutorials/alphagenome-upstream/colabs/`)
- Got AlphaGenome API key at https://alphagenome.google/api (free, non-commercial)
- Wrote biology primer for future-me

### What I learned
- AlphaGenome is a Python gRPC client wrapping a remote API, not a local model
- The `dna_model.create()` factory is the entry point
- Three core operations: predict (track values), score_variant (ref vs alt delta), ISM (in-silico mutagenesis)
- Free tier allows ~1000s of predictions, not millions (use Atlas for genome-wide queries)

### What confused me
- ~~v1 vs v2 API confusion~~ resolved: SDK is current, tutorials use current classes (`genome.Interval`, `TrackData`)
- Sequence length must be 4,096 ≤ L ≤ 1,048,576 — anything shorter needs padding (the model centers on a window)

### What I learned reading the source
- **Two clients exist**: `dna_client` (live API) and `atlas` (pre-computed scores). Atlas is the right tool for genome-wide scans.
- **Hard limits**: `MAX_VARIANT_SCORERS_PER_REQUEST = 20`, `MAX_ISM_INTERVAL_WIDTH = 10` (per chunk, internal batching), 4 sequence lengths (16K/100K/500K/1M).
- **18 recommended variant scorers** — `_ACTIVE` variants are tuned for variant scoring and often more sensitive than base versions.
- **5 model folds** (`ModelVersion.FOLD_0` through `FOLD_3` plus `ALL_FOLDS`). Default uses ensemble; benchmarking wants individual folds.
- **Organism**: human + mouse, but mouse is limited (most scorers don't have mouse tracks).
- **Atlas chunks intervals into 32 bp pieces** and runs them in parallel with 10 workers.

### What I built (in addition to README + primer)
- `scripts/smoke_test.py` — verify API key in ~30s
- `scripts/atlas_smoke_test.py` — verify Atlas access, see what data comes back
- `scripts/quickstart_script.py` — full quick_start.ipynb flow as a script
- `scripts/ism_example.py` — ISM walkthrough with heatmap
- `data/gencode.v46.annotation.gtf.gz.feather` (318 MB, pre-fetched)
- `notes/sdk-cheatsheet.md` — distilled SDK knowledge
- `figures/` directory for outputs

3. Pipeline crashed mid-merge (OOM during in-memory accumulation of AnnData chunks)
4. Atlas SDK bug confirmed: all approaches fail at ~47% with `gRPC semaphore leaked` warning
5. **Pivot to live API (`dna_client.score_variant`)** — bypass Atlas SDK bug
6. Wrote `scripts/benchmark_scn1a_live_api.py` — score 591 ClinVar variants (216 pathogenic splicing + 375 benign intronic)
7. **Result: AUPRC = 0.9833 for SPLICE_SITES, top-5% precision = 100%, AUROC = 0.9896**
8. All 591 API calls succeeded (no failures), runtime ~6.5 minutes
9. Bimodal score distribution: pathogenic ~1.0, benign ~0.05, almost no overlap

### Status
- **Phase 1-3: DONE**
- Methods paper is realistic, not aspirational
- Need: paper draft, SpliceAI comparison, novel candidate identification, outreach

### Next actions
1. Run SpliceAI on same variants for comparison (if accessible)
2. Identify novel candidates: ClinVar VUS in SCN1A that score high
3. Begin paper draft: intro (Dravet background), methods (pipeline), results (this benchmark), discussion
4. Outreach to Carvill / Sparber labs

## Day 2

### What I did
- Picked Dravet syndrome / SCN1A as the target disease
- Wrote `PROJECT_PLAN.md` with phased delivery and 90-day timeline
- Curated `data/known_variants.py` with the Carvill 2018 / Sparber 2023 known pathogenic variants
- Downloaded ClinVar VCF (185 MB) + tabix index
- Wrote `scripts/clinvar_scn1a_local.py` — pulls all SCN1A ClinVar variants, labels them pathogenic / splicing-related
- Result: 5,276 SCN1A records, 216 pathogenic + splicing-related variants as gold-standard positive controls

### What I learned
- ClinVar VCF uses contig names without "chr" prefix (`2` not `chr2`)
- ClinVar's GENEINFO field uses **pipe** separator, not comma: `SCN1A:6323|LOC102724058:102724058`
- ClinVar's INFO values can contain commas inside (e.g., CLNDISDB has comma-separated database IDs), so naive `split(",")` breaks things
- pysam.TabixFile.fetch() returns an **iterator** that gets consumed by `len(list(...))` — always materialize to a list first

### Bugs I hit and fixed
- v1 ClinVar fetch via Entrez returned only 21 records (using wrong query syntax)
- v2 ClinVar REST API returned 404 (deprecated)
- Tried remote tabix query (HTTPS) — not supported by pysam
- Fixed: download VCF + index locally, query by region
- Fixed: comma vs pipe separator bug in GENEINFO parsing
- Fixed: tabix iterator consumption bug

### Next actions
1. Build `scripts/atlas_scn1a.py` — query AlphaGenome Atlas for the SCN1A locus
2. Cross-reference ClinVar pathogenic variants with Atlas scores
3. Build benchmark: AUPRC per splicing modality
4. Identify novel candidate regulatory variants
5. Begin writing the paper draft
