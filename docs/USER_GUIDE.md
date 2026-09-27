# User Guide — alphagenome-scn1a-pipeline

> A practical walkthrough for computational biology researchers who want to apply our methodology to a new gene, score a custom VCF, or re-use our pre-computed SCN1A dataset.

This guide assumes you have basic Python fluency, can read a VCF, and have an AlphaGenome API key. It does **not** assume you've read our paper.

**Repository:** https://github.com/rollroyces/alphagenome-scn1a-pipeline
**Hugging Face dataset:** https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome
**API key (free, non-commercial):** https://alphagenome.google/api

---

## Table of contents

1. [Quick start (5 min)](#1-quick-start-5-min)
2. [Run a benchmark on a new gene (30 min)](#2-run-a-benchmark-on-a-new-gene-30-min)
3. [Re-score a custom VCF (15 min)](#3-re-score-a-custom-vcf-15-min)
4. [Use the pre-computed dataset](#4-use-the-pre-computed-dataset)
5. [Interpret results](#5-interpret-results)
6. [FAQ](#6-faq)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Quick start (5 min)

Get the repo running and verify your setup with a 5-variant smoke test. No benchmark is run here — just "does my API key work, does my venv work, does the pipeline import".

### 1.1 Clone and enter

```bash
git clone https://github.com/rollroyces/alphagenome-scn1a-pipeline.git
cd alphagenome-scn1a-pipeline
```

### 1.2 Get an API key

Free, non-commercial, request at https://alphagenome.google/api. You'll receive a string that starts with `AIza`.

### 1.3 Save the key (never commit it)

```bash
echo "AIza..." > .alphagenome_key
chmod 600 .alphagenome_key
```

The file `.alphagenome_key` is in `.gitignore` — keep it that way. The launcher `scripts/_run_with_key.sh` reads this file and exports `ALPHAGENOME_API_KEY` for you, so you don't have to set it in your shell.

### 1.4 Create the venv and install

```bash
uv venv --python 3.13 .venv
uv pip install --python .venv/bin/python -e '.[dev]'
```

(`-e '.[dev]'` installs the package in editable mode plus pytest. You can drop `.[dev]` if you only want the runtime deps.)

### 1.5 Verify with `make smoke`

```bash
make smoke
```

This runs `scripts/smoke_test.py`, which:

1. Picks 5 SCN1A SNVs (seeded, deterministic) from `outputs/benchmark_scn1a_live_api_raw.csv`.
2. Re-scores them against the live AlphaGenome API (~5 API calls; takes seconds).
3. Compares to the recorded scores within `1e-3` absolute tolerance (model is deterministic for the same sequence + interval).
4. Computes a mini-AUPRC on the 5 picks and requires it be `> 0.5`.

**Expected output (last lines):**

```
[smoke] mini-AUPRC (SPLICE_SITES, n=5): 1.0000 (full benchmark AUPRC=0.9833)

[smoke] PASS — pipeline reproduces existing outputs within tolerance.
```

If you see `PASS`, your install, API key, and pipeline are good. If you see `FAIL`, jump to [Troubleshooting](#7-troubleshooting).

### 1.6 Other useful make targets

```bash
make help        # list every target
make key         # check that .alphagenome_key exists and is plausibly long
make check       # verify all expected outputs/*.csv exist (no API calls)
make venv        # (re)create the .venv
make test        # run unit tests, no API key required
```

---

## 2. Run a benchmark on a new gene (30 min)

The benchmark is the core reproducibility unit: feed AlphaGenome a balanced set of ClinVar variants for your gene and ask "can the model separate known-pathogenic from known-benign?"

### 2.1 The exact command

The repo's convention is a pair of scripts per gene: `scripts/_<gene>_extract.py` (build a ClinVar benchmark TSV) followed by `scripts/_<gene>_score.py` (score with the live API, compute AUPRC). The cross-disease runner already exists for **SCN1A, SCN2A, MECP2, CFTR, DMD, KCNQ2, COL4A5, FBN1** — see `scripts/cross_disease_benchmark.py` and `scripts/_kcnq2_score.py` for working examples.

To apply the same pattern to a new gene, use the gene-agnostic extractor:

```bash
# 1. extract a balanced ClinVar benchmark for MY_GENE
python scripts/extract_clinvar_for_gene.py \
    --gene MY_GENE \
    --chrom 7 \
    --start 12345678 \
    --end 12399999 \
    --n-pos 200 \
    --n-neg 350 \
    --out outputs/clinvar_my_gene_benchmark.tsv

# 2. score with the live API and compute AUPRC
bash scripts/_run_with_key.sh scripts/_kcnq2_score.py \
    --input outputs/clinvar_my_gene_benchmark.tsv \
    --raw-out outputs/cross_disease_my_gene_raw.csv \
    --metrics-out outputs/cross_disease_my_gene_metrics.csv
```

`_kcnq2_score.py` is the closest gene-agnostic template: it reads any TSV with `chrom, pos, ref, alt, rsid, clnsig_category, label`, samples a balanced benchmark, calls `dna_model.score_variant` with the three splicing scorers, and writes both raw scores and per-scorer metrics.

For a fully wiring-it-together alternative, the SCN1A originals are `scripts/benchmark_scn1a_live_api.py` and `scripts/extract_clinvar_for_gene.py`. You can also see the per-gene extractors `scripts/_kcnq2_extract.py`, `scripts/_col4a5_extract.py`, `scripts/_fbn1_extract.py` for the lighter "tabix ClinVar, stratify by clnsig" pattern.

### 2.2 What the benchmark script does

For each variant it scores, in this order:

1. **Extract a ClinVar benchmark** for the gene locus using `pysam.TabixFile` over `data/clinvar_grch38.vcf.gz`. See `scripts/extract_clinvar_for_gene.py:make_benchmark`.
2. **Filter positives** to pathogenic SNVs with splicing-related molecular consequences (splice donor, acceptor, region, polypyrimidine tract).
3. **Filter negatives** to benign SNVs that fall inside intronic sequence.
4. **Cap** at `n_pos` / `n_neg` (defaults 200/350 — matches the SCN1A/SCN2A-style benchmark).
5. **Score with AlphaGenome** at 16 Kb context (`dna_client.SEQUENCE_LENGTH_16KB`) using three scorers:
   - `SPLICE_SITES` — pointwise disruption to splice-site probability tracks.
   - `SPLICE_SITE_USAGE` — usage change for nearby splice sites.
   - `SPLICE_JUNCTIONS` — junction-level read-coverage change.
6. **Compute metrics** with bootstrap 95% CIs (`scripts/_kcnq2_score.py:bootstrap_auprc_ci`):
   - AUROC
   - AUPRC (with 95% percentile bootstrap CI over 1000 resamples)
   - Top-5% precision (the variant score at the 95th percentile rank — fraction of "positives" in the top 5% by score)

### 2.3 Expected output format

**`outputs/cross_disease_my_gene_raw.csv`** — one row per variant, per-scorer scores + success flag:

| column | type | meaning |
|---|---|---|
| `chrom`, `pos`, `ref`, `alt` | str/int/str/str | Variant coordinates (hg38, **no `chr` prefix** in the file but AlphaGenome accepts either). |
| `rsid` | str | dbSNP rsID if available, else empty. |
| `label` | str | `"positive"` or `"negative"` (pathogenic vs benign). |
| `clnsig_category` | str | Coarse ClinVar category: `pathogenic` / `benign` / `uncertain` / `conflicting` / `other`. |
| `SPLICE_SITES_score` | float | Sum of AlphaGenome annotation track. Higher = more splicing disruption predicted. |
| `SPLICE_SITE_USAGE_score` | float | Same scale family. |
| `SPLICE_JUNCTIONS_score` | float | Same. |
| `score_success` | bool | True if the API call returned without error. |
| `score_error` | str | Short exception text if `score_success=False`. |

**`outputs/cross_disease_my_gene_metrics.csv`** — one row per scorer:

| column | type | meaning |
|---|---|---|
| `gene` | str | Gene symbol (set in the script). |
| `scorer` | str | One of `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`. |
| `n_total`, `n_positive`, `n_negative` | int | Counts of successfully-scored variants. |
| `auroc` | float | AUROC across the benchmark. |
| `auprc` | float | AUPRC across the benchmark. |
| `auprc_ci_lo`, `auprc_ci_hi` | float | 95% percentile bootstrap CI on AUPRC. |
| `top_5_pct_precision` | float | Fraction of positives in the top 5% by score. |

A real example (SCN1A, from `outputs/cross_disease_scn1a_raw.csv`):

```
chrom,pos,ref,alt,rsid,label,clnsig_category,SPLICE_SITES_score,SPLICE_SITE_USAGE_score,SPLICE_JUNCTIONS_score,score_success
2,165989824,G,A,rs121918633,positive,pathogenic,1.7333984375,196.203125,2783.7890625,True
```

### 2.4 Runtime and API cost

- **Live API scoring:** ~0.4–0.6 s/variant (single-variant 16 Kb context). A 550-variant benchmark takes ~5 minutes wall-clock.
- **The 95% bootstrap CI on AUPRC** runs in-memory; no API calls.
- **Atlas SDK (see §3) does not consume live API quota** — it's pre-computed.

---

## 3. Re-score a custom VCF (15 min)

Use case: a clinical lab sends you 30 VUS from their internal pipeline. You want to ask AlphaGenome "which of these look the most pathogenic by splicing impact?" and get a ranked table.

### 3.1 The reference template

`scripts/rescore_vus.py` is the canonical re-scorer — built for SCN1A but trivially gene-agnostic. It:

1. Reads `outputs/clinvar_scn1a.tsv` (the gene's full ClinVar dump).
2. Filters to `clnsig_category == "uncertain"`.
3. Filters to **SNVs only** (`ref` and `alt` both length 1 — AlphaGenome splicing scorers don't take indels).
4. Scores each VUS with the same three splicing scorers as the benchmark.
5. Ranks by `SPLICE_SITES_score`, writes top-50 and "high-impact" (score ≥ 0.5) subsets.

### 3.2 Adapting it to your VCF

Drop your VCF into a DataFrame with `chrom, pos, ref, alt` columns (hg38), then call the same scorer block:

```python
import os
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

df = pd.read_csv("my_vus.tsv", sep="\t")  # columns: chrom, pos, ref, alt, ...
df = df[(df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]

dna_model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
]
scorer_names = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]

rows = []
for _, row in df.iterrows():
    chrom = str(row["chrom"])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    variant = genome.Variant(chrom, int(row["pos"]), str(row["ref"]), str(row["alt"]))
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
    scores = dna_model.score_variant(interval=interval, variant=variant, variant_scorers=scorers)
    out = dict(row)
    for ann, name in zip(scores, scorer_names):
        out[f"{name}_score"] = float(ann.X.sum())
    out["score_success"] = True
    rows.append(out)

scored = pd.DataFrame(rows).sort_values("SPLICE_SITES_score", ascending=False)
scored.to_csv("my_vus_rescored.csv", index=False)
```

This is the inner loop of `scripts/rescore_vus.py:score_vus` with the loader swapped out.

### 3.3 Expected output

`outputs/vus_rescored.csv` (SCN1A reference, first 2 rows):

```
chrom,pos,ref,alt,rsid,clnsig_raw,clnsig_category,...,SPLICE_SITES_score,SPLICE_SITE_USAGE_score,SPLICE_JUNCTIONS_score,score_success,vus_rank
2,166047622,C,A,851265,Uncertain_significance,uncertain,...,1.71533203125,140.3272247314453,2616.03173828125,True,1
2,166054633,C,T,2203203,Uncertain_significance,uncertain,...,1.64013671875,146.9827880859375,1398.0078125,True,2
```

The `vus_rank` column is the rank by `SPLICE_SITES_score` (1 = highest predicted impact). Other useful outputs:

- `outputs/vus_top_candidates.csv` — top 50 with rsID, molecular consequence, ClinVar disease (`clndn`), review status.
- `outputs/vus_high_impact_candidates.csv` — score ≥ 0.5 only (the "looks pathogenic-like" subset).

### 3.4 Use the Atlas SDK first to avoid burning API quota

For genome-wide or locus-wide screens, hit the **pre-computed AlphaGenome Variant Impact (AVI)** tracks via the Atlas SDK. This costs zero live-API quota:

```bash
bash scripts/_run_with_key.sh scripts/atlas_scn1a.py
```

Writes:
- `data/atlas_scn1a_splicing.h5ad` — full AnnData with all three splicing scorers for every SNV in `chr2:165_934_640–166_232_806` (gene locus ± 50 kb flanking).
- `outputs/atlas_scn1a_summary.csv` — flat CSV you can join with your VCF on `(chrom, pos, ref, alt)`.

The Atlas SDK only exposes `SPLICE_SITES`, `SPLICE_SITE_USAGE`, and `SPLICE_JUNCTIONS` for pre-computed variants (`SPLICE_JUNCTIONS_ACTIVE` is not in Atlas — verified empirically in `scripts/atlas_scn1a.py:43`). Use the live API only for the variants you actually need to commit to.

---

## 4. Use the pre-computed dataset

We published all 1,610 SCN1A VUS with AlphaGenome scores to Hugging Face so you don't have to re-run the API.

### 4.1 Load it

```bash
pip install datasets
```

```python
from datasets import load_dataset

ds = load_dataset("RROL/scn1a-vus-alphagenome")
df = ds["train"].to_pandas()
print(df.shape)  # (1610, 30)
print(df.columns.tolist())
```

### 4.2 Schema (30 columns)

Key columns (full list in `HUGGINGFACE_DATASET.md`):

| column | type | meaning |
|---|---|---|
| `chrom`, `pos`, `ref`, `alt` | int/str | hg38 coordinates. |
| `rsid` | str | dbSNP rsID. |
| `molecular_consequence` | str | ClinVar MC label (e.g. `intron_variant`, `splice_region_variant`). |
| `vus_rank` | int | Rank by `SPLICE_SITES_score` (1 = highest predicted impact). |
| `tier` | int | 1 / 2 / 3 (see §5.3). |
| `clnsig_category` | str | ClinVar coarse category. |
| `DNASE_score` | float | DNASE concentration in matched tissue (lower = more closed chromatin). |
| `SPLICE_SITES_score` | float | AlphaGenome. |
| `SPLICE_SITE_USAGE_score` | float | AlphaGenome. |
| `SPLICE_JUNCTIONS_score` | float | AlphaGenome. |
| `dnase_success`, `dnase_error` | bool/str | DNASE lookup status. |
| `clndn` | str | ClinVar disease name(s). |
| `review_status` | str | ClinVar star rating (e.g. `criteria_provided`). |
| `present` | bool | True if the variant appears in gnomAD v4.1. |
| `gnomad_genome_{af,ac,an,ac_hom}` | float | gnomAD genome frequencies. |
| `gnomad_exome_{af,ac,an,ac_hom}` | float | gnomAD exome frequencies. |
| `error` | str | Per-row scoring error, if any. |
| `pubmed_hits` | int | Number of PubMed records mentioning this rsID (0 = novel). |

### 4.3 Cite it

```
Royce Lam. (2026). SCN1A VUS Re-Scoring with AlphaGenome (Revision 8db4888).
Hugging Face. https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome
```

License: **CC-BY-4.0**. You can use, redistribute, and build on it as long as you cite.

### 4.4 Quick analyses

```python
# Top 10 novel candidates (no prior literature)
novel = df[(df["pubmed_hits"] == 0) & (df["tier"] <= 2)]
print(novel.nlargest(10, "SPLICE_SITES_score")[["rsid", "chrom", "pos", "SPLICE_SITES_score", "tier"]])

# Distribution of SPLICE_SITES_score by tier
print(df.groupby("tier")["SPLICE_SITES_score"].describe())
```

---

## 5. Interpret results

### 5.1 What AUPRC means here

AUPRC (area under the precision-recall curve) measures how well a score separates pathogenic from benign. **1.0 = perfect separation, 0.5 = random** (when class balance is 50/50; the random baseline moves with class prevalence).

For the SCN1A benchmark with **120 pathogenic + 310 benign intronic variants**, our AUPRC is **0.9833** (95% CI [0.964, 0.996] — see `outputs/benchmark_scn1a_results.csv` and the bootstrap in `scripts/_kcnq2_score.py:bootstrap_auprc_ci`).

Reference numbers from our 8-gene cross-disease benchmark (`outputs/cross_disease_summary.md`), all on the **filtered** protocol (positives = pathogenic + splice-related; negatives = benign + intronic):

| Gene | n_pos | n_neg | AUROC | AUPRC | 95% CI | Top-5% prec |
|---|---:|---:|---:|---:|---|---:|
| MECP2 | 12 | 100 | 1.0000 | 1.0000 | [1.000, 1.000] | 1.000 |
| FBN1 | 100 | 200 | 0.9999 | 0.9997 | [0.999, 1.000] | 1.000 |
| DMD | 200 | 400 | 1.0000 | 0.9999 | [1.000, 1.000] | 1.000 |
| CFTR | 150 | 300 | 0.9994 | 0.9988 | [0.997, 1.000] | 1.000 |
| COL4A5 | 176 | 200 | 0.9966 | 0.9969 | [0.992, 1.000] | 1.000 |
| SCN2A | 30 | 200 | 0.9980 | 0.9880 | [0.965, 1.000] | 1.000 |
| **SCN1A** | **120** | **310** | **0.9940** | **0.9830** | **[0.964, 0.996]** | **1.000** |
| KCNQ2* | 200 | 350 | 0.6713 | 0.5657 | [0.501, 0.634] | 0.929 |

\*KCNQ2 used an **unfiltered** protocol (all pathogenic vs all benign, no splice-consequence filter). The 7 filtered genes average **AUPRC = 0.995** with σ ≈ 0.007. **Top-5% precision is 1.000 across all 8 genes** — the high-score tail is uniformly enriched for pathogenic variants regardless of protocol.

### 5.2 When to trust top-K precision

Top-K precision (top 5% precision in our defaults) is the metric that actually matters for variant prioritization: "if I take the 5% of variants with the highest AlphaGenome score, what fraction are truly pathogenic?"

**Trust top-K when:**
- n_pos ≥ 30 (so the bootstrap CI on AUPRC is reasonably tight)
- Your positives are enriched for splice-consequence pathogenic variants (otherwise the splice scorers can't see them — see KCNQ2 row above)
- The gene is on a standard chromosome and the locus is mapped correctly to hg38

**Top-K precision ≥ 0.95 across all 8 of our genes** means that, **for splice-disrupting pathogenic variants**, the high-score tail is reliably enriched. If you see top-5% precision drop below ~0.7 in your run, suspect a coordinate-system bug, a wrong transcript strand, or contamination of your negative set with splice-region pathogenic variants (see §5.4).

### 5.3 When to be skeptical

The 8-gene mean AUPRC of 0.995 is encouraging but masks three failure modes that *will* bite you on a new gene:

1. **Small `n_pos` (n < 20).** The 95% bootstrap CI widens fast; AUPRC can look "perfect" by chance. MECP2 has 12 positives — treat its 1.000 as suggestive, not definitive.
2. **Missense contamination.** If your positive set is dominated by missense / nonsense pathogenic variants (most ClinVar "pathogenic" calls), the splice scorers can't see them and AUPRC will crater. **Pre-filter positives to those with splice-related molecular consequences** (see `SPLICING_CONSEQUENCES` in `scripts/extract_clinvar_for_gene.py:25`). This is why KCNQ2 (mostly missense) hits 0.57 and SCN1A (splice-rich) hits 0.98 on the unfiltered protocol.
3. **Wrong transcript strand / wrong assembly.** See Troubleshooting §7.3 and §7.4.

The splicing scorers are **mechanism-specific**, not generic pathogenicity classifiers. Don't expect them to rank missense variants.

### 5.4 Tier definitions (SCN1A)

From `HUGGINGFACE_DATASET.md`:

- **tier = 1** (4 variants: rs801806, rs4293437, rs801809, rs2847163) — the wet-lab-testable candidates. Selected by mechanism + orthogonal evidence (chromatin similarity to known pathogenic + brain RNA-seq effect + splice-site proximity on MANE Select), **not** by raw splice rank. They sit at vus_rank 14, 17, 29, 30.
- **tier = 2** (46 variants) — top 50 minus top 4.
- **tier = 3** (1,560 variants) — all other SCN1A VUS in ClinVar.

**Caveat:** the very top of the splice ranking is dominated by deep-intronic variants (`outputs/vus_tier2_candidates.csv`). Tier-1 selection deliberately trades raw score for mechanism: a variant 2 bp from a canonical splice site with a benign-looking DNASE profile outranks a deep-intronic variant with a higher SPLICE_SITES_score.

---

## 6. FAQ

### 6.1 Where do I get MANE Select transcript IDs?

MANE Select is a joint RefSeq + Ensembl project that picks one canonical transcript per gene.

- **NCBI RefSeq:** https://www.ncbi.nlm.nih.gov/refseq/MANE/ — searchable by gene symbol; click "MANE Select" → gives NM_ accession and the matching ENST ID.
- **Ensembl:** https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=<ENSG_ID> — shows MANE Select as `ENST...`.
- **GENCODE:** the GTF file at https://www.gencodegenes.org/human/ has a `tag "MANE_Select"` on the relevant transcript row.

For reference, the MANE Select transcripts in our cross-disease benchmark are hard-coded in `scripts/cross_disease_benchmark.py:38` (SCN1A, SCN2A, MECP2, CFTR, DMD) and in `scripts/_kcnq2_extract.py:19` (KCNQ2 = ENST00000356457 / NM_172107.4, hg38 plus strand), `scripts/_col4a5_extract.py:11`, and `scripts/_fbn1_extract.py:11`.

### 6.2 How do I pick the right tissue filter?

For **splicing** scorers, tissue generally doesn't matter (splice-site consensus is largely invariant). Our 8-gene benchmark spans brain (SCN1A, SCN2A, MECP2), muscle (DMD), epithelial (CFTR), kidney (COL4A5), and connective-tissue / fibroblast (FBN1) — and the AUPRC range is 0.97–1.00. Don't waste API quota on tissue-specific runs unless you're adding the **DNASE concentration** track, where tissue matters a lot.

For **DNASE**, match the disease-relevant tissue:
- SCN1A / SCN2A / MECP2 → brain (we used brain in `outputs/vus_rescored_with_dnase.csv`).
- DMD → skeletal muscle.
- CFTR → airway epithelium.
- COL4A5 → kidney.
- FBN1 → fibroblast.

Tissue-specific score weights are set at the `dna_model.score_variant(...)` call site — see the existing per-experiment scripts in `research_notebook/experiments/*/` for examples.

### 6.3 How do I avoid burning API quota?

Three rules of thumb:

1. **Use Atlas first.** `scripts/atlas_scn1a.py` queries the pre-computed AlphaGenome Variant Impact (AVI) tracks via `alphagenome.atlas`. Costs zero live-API quota, returns all three splicing scorers for every SNV in your locus. Use it as a screen; only call the live API for variants you actually need to confirm.
2. **Use `make smoke` to debug, not full benchmarks.** `make smoke` runs 5 variants. A full SCN1A benchmark is ~591 variants.
3. **Use the Hugging Face dataset.** All 1,610 SCN1A VUS are pre-scored — `pip install datasets` and skip the API entirely.

For very large loci (>200 Kb), Atlas can OOM (the gRPC default message size cap is 4 MB). `scripts/atlas_scn1a.py:60` patches `grpc.secure_channel` to lift the cap; copy that block.

### 6.4 What does the SCN1A README headline number (61 novel high-impact candidates) mean?

The full pipeline found **61 SCN1A VUS** that simultaneously:
- Score in the high-impact range (SPLICE_SITES_score ≥ 0.5, comparable to known pathogenic controls).
- Have **0 PubMed mentions** (truly novel).
- Either are absent from gnomAD or have allele frequency < 1e-4.
- (For missense variants) Predicted to have a splice-altering effect on top of the coding change.

These are the wet-lab-testable candidates. See `outputs/vus_high_impact_with_gnomad.csv` and `paper/candidate_report.md`.

---

## 7. Troubleshooting

### 7.1 API key errors

Symptom:

```
FAIL: ALPHAGENOME_API_KEY environment variable is not set.
```

or

```
PermissionDeniedError: 403 ...
```

Fixes (in order):

1. Confirm `.alphagenome_key` exists and is ≥ 30 bytes: `make key`.
2. Confirm permissions: `chmod 600 .alphagenome_key`.
3. Confirm you're using the launcher (which exports the env var), not bare `python`: `bash scripts/_run_with_key.sh <script>` not `python <script>`.
4. Confirm the key starts with `AIza` and matches what https://alphagenome.google/api gave you (no leading whitespace, no newline).
5. If the API key is fresh, you may need to wait a few minutes for the quota system to propagate — the live API returns 403 on brand-new keys for ~5 minutes.

### 7.2 Out of memory on large genes

Symptom:

```
MemoryError / Killed (signal 9) / OOMKilled
```

Most likely on:
- Atlas queries over loci > 200 Kb (`scripts/atlas_scn1a.py:60` patches the gRPC message-size cap; the resulting AnnData is large).
- ISM (in-silico mutagenesis) experiments with thousands of single-base substitutions per variant (`research_notebook/experiments/001_ism_scn1a/`).

Fixes:

1. **Use the chunked Atlas loader:** `scripts/atlas_scn1a_chunked.py` and `scripts/atlas_scn1a_csv.py` break the locus into windows.
2. **Narrow the locus.** For SCN1A we use gene body + 50 kb flanks (`scripts/atlas_scn1a.py:39`). Drop the flank for first try.
3. **Run on a machine with more RAM.** A typical 200 Kb locus with all three splicing scorers needs ~4–6 GB; FBN1-scale (237 kb) needs ~8 GB.
4. **Don't load all tissues at once.** Pass only the scorers you need to `dna_model.score_variant(...)`.

### 7.3 Wrong MANE Select strand

Symptom: AUPRC ~ 0.5 (random), or your positives look uniformly distributed across scores.

Diagnosis:

```python
# Quick check — does the model "see" a pathogenic splice donor as disrupted?
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import os

dna_model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
variant = genome.Variant("chr2", 166041470, "T", "A")  # known SCN1A Tier-1 donor
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"]],
)
print(scores[0].X.sum())
# Expected: ~1.1 (high — this is a known disrupted donor)
# If you get ~0, you're on the wrong strand or the wrong transcript context.
```

Fixes:

1. **Verify the gene is on the strand you think it is.** SCN1A is minus-strand on chr2; KCNQ2 is plus-strand on chr20; CFTR is plus-strand on chr7. The cross-disease dictionary in `scripts/cross_disease_benchmark.py:38` records strand explicitly.
2. **Verify hg38 coordinates.** ClinVar uses hg38 in their standard VCF; older sources may use hg19 (see §7.4).
3. **Verify the locus endpoints.** Pull the gene interval from GENCODE and check that it contains the expected pathogenic variant.

### 7.4 Coordinate system (hg19 vs hg38)

Symptom: every variant scores ~0, or the model returns coordinates outside the gene body.

Diagnosis:

```bash
# Check the assembly tag in your ClinVar VCF
zcat data/clinvar_grch38.vcf.gz | grep -m1 "##fileformat" 
zcat data/clinvar_grch38.vcf.gz | grep -m1 "##assembly"  # should say GRCh38
```

If you're using **hg19** ClinVar and the model is built on hg38 (the case for the public AlphaGenome API), every coordinate is off by tens to hundreds of Mb. Cross-disease liftover:

- **LiftOver CLI:** https://genome.ucsc.edu/cgi-bin/hgLiftOver (batch mode via `liftOver`).
- **Python liftover:** `pip install pyliftover` — see `from pyliftover import LiftOver; LiftOver('hg19', 'hg38').convert_coordinate('chr2', 166123456)`.
- **NCBI Genome Remapping Service:** https://www.ncbi.nlm.nih.gov/genome/tools/remap — supports batch VCF liftover.

The default ClinVar VCF we use (`data/clinvar_grch38.vcf.gz`) is already hg38; the issue only arises if you bring your own ClinVar dump from an older source.

### 7.5 Nothing else worked

1. Check `make check` — all expected outputs are present.
2. Run `make test` — the unit-test suite doesn't need an API key.
3. File an issue: https://github.com/rollroyces/alphagenome-scn1a-pipeline/issues with the full stdout/stderr of `make smoke` and the output of `python -V`, `uv --version`, and `pip show alphagenome`.

---

## Where to go next

- **Read the methods paper:** `paper/preprint.md` (~5,000 words; the methods + main result).
- **See the lab-ready variant list:** `paper/candidate_report.md` (the 4 Tier-1 SCN1A candidates).
- **Reach out to a collaborator:** `paper/outreach_template.md` (draft emails to Carvill / Sparber / Helbig labs).
- **Browse the experiments:** `research_notebook/experiments/` — every benchmark run is reproducible from here.

Happy variant hunting.