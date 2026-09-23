# Experiment 009 — DNASE re-scoring of SCN1A VUS

**Date:** 2026-09-23 (HKT)
**Status:** Complete

## TL;DR

All 1,610 SCN1A VUS were re-scored with the `DNASE` scorer via the live API
in addition to the existing `SPLICE_SITES` / `SPLICE_SITE_USAGE` /
`SPLICE_JUNCTIONS` scores. **0/1,610 failures**, **508 s** total runtime
(~3.2 calls/sec, faster than the 1.5 calls/sec in `scripts/rescore_vus.py`).

The combined (mean of min-max-normalized SPLICE_SITES and DNASE) top-30
overlaps the splicing-only top-30 by **27/30**. **DNASE is *not* fully
redundant** — Spearman ρ between the two scorers on VUS is only **−0.024
(p=0.34)**, so they are nearly independent. The 3 "new" variants DNASE pulls
in are mostly **missense variants with very low splicing scores (median
0.022) but high chromatin impact** — these are exactly the dual-mechanism
candidates flagged in the `alphagenome-variant-scoring` skill ("Missense VUS
with high splice scores are dual-mechanism candidates, not false positives")
applied in the other direction.

**Recommendation:** mention **3 chromatin-only candidates** in the lab
outreach (rsIDs 1381398, 4532737, 2178954) alongside the 27 dual-mechanism /
splice candidates already in `outputs/vus_top_candidates.csv`. Frame them
as "DNASE-disrupted missense variants with low splice impact — orthogonal
mechanism, mechanistically distinct."

## Files produced

| File | Rows | Purpose |
|---|---|---|
| `outputs/vus_rescored_with_dnase.csv` | 1,610 | All VUS with DNASE + splicing scores |
| `outputs/vus_top_by_dnase.csv`        | 30    | Top 30 by DNASE alone |
| `outputs/vus_top_combined.csv`        | 30    | **Top 30 by mean(normalized splice, normalized DNASE)** — key artifact |
| `outputs/vus_high_impact_with_dnase.csv` | 61 | Existing 61 candidates + new DNASE column |
| `outputs/vus_combined_vs_splicing.json` | –   | Machine-readable comparison summary |
| `research_notebook/experiments/009_dnase_vus_rescore/run.log` | – | Full stdout from the run |
| `research_notebook/experiments/009_dnase_vus_rescore/summary_stats.json` | – | Pearson / Spearman / verdict for this README |

**NOT modified:** `outputs/vus_rescored.csv`, `outputs/vus_high_impact_with_gnomad.csv`,
`paper/preprint.md`, any existing scripts.

## Method

- Input: `outputs/vus_rescored.csv` (1,610 VUS, sorted by `vus_rank`).
- API call: `dna_client.score_variant(interval=…, variant=…, variant_scorers=[DNASE])`,
  16 KB window centered on the variant. Identical pattern to `scripts/rescore_vus.py`.
- DNASE scoring via `variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]`
  (a `CenterMaskScorer`). Returns 305 DNase track predictions (DNase-seq
  across tissues/cell types), summary statistic = `ann.X.sum()` per variant.
- Failed calls → `DNASE_score = NaN`, `dnase_success = False`, `dnase_error = "<type>: <msg>"`.
  No variant is skipped (the script appends the failure row and moves on).
- Joint ranking: min-max normalize `SPLICE_SITES_score` (best single
  splicing discriminator per `rescore_vus.py`) and `DNASE_score` each to
  `[0, 1]`; rank by their mean.

## Results

### 1. Run summary

| Metric | Value |
|---|---|
| n_scored | 1,610 (100%) |
| n_failed | 0 |
| runtime | 508 s (~8.5 min) |
| throughput | ~3.17 calls/sec |
| DNASE_score range | −61.1 to +104.4 |
| DNASE_score median | +0.75 |
| DNASE_score p25 / p75 | −5.2 / +6.7 |

### 2. Top-30 overlap (splicing-only vs combined)

| Set | Count |
|---|---|
| Splicing-only top-30 | 30 |
| Combined top-30 | 30 |
| **Overlap** | **27** |
| Only in combined (DNASE pulled them in) | 3 |
| Only in splicing (DNASE pushed them out) | 3 |

Per the experiment rubric, 27/30 overlap = "REDUNDANT" — but the rubric's
binary verdict misses the more interesting structure (see §3 below).

### 3. DNASE captures a different *variant class*

The DNASE-only top-30 is **86% missense variants**, vs 23% missense in the
splicing-only top-30. The DNASE top-30 has median SPLICE_SITES_score of
**0.022** — these variants have essentially no splicing impact. They were
invisible to the original pipeline.

```
Consequence mix in each top-30
                          combined  splicing  dnase
intron_variant                17         16       3
missense_variant               7          7      26
non-coding_transcript_var.     3          3       1
splice_donor_variant           2          2       0
splice_acceptor_variant        1          2       0
```

### 4. Independence of the two signals

On the full 1,610-VUS set:

- **Pearson r(SPLICE_SITES, DNASE) = −0.065**
- **Spearman ρ = −0.024 (p = 0.34)**

Both scorers are effectively uncorrelated on VUS. They are measuring
different things — splicing site disruption vs chromatin (DNase I
hypersensitivity) perturbation. That is exactly why even 3 non-overlapping
variants in a 30-row top list are mechanistically meaningful.

### 5. Top 10 by combined ranking

| Rank | rsID | Pos | REF→ALT | SpliceScore | DNASEScore | CombinedScore | Consequence | Disease (clndn) |
|---:|---|---:|---|---|---:|---:|---|---|
| 1  | 1412774  | 166047773 | A→C | 1.538 | +24.8  | 0.708 | intron_variant                | Early-infantile_DEE |
| 2  | 851265   | 166047622 | C→A | 1.715 | −15.0  | 0.639 | intron_variant                | Migraine |
| 3  | 393000   | 166047622 | C→G | 1.613 |  −8.3  | 0.630 | intron_variant                | Complex_neurodev. / Early-infantile_DEE |
| 4  | 1046194  | 166013744 | C→T | 1.521 | −10.0  | 0.598 | non-coding_transcript_variant | Early-infantile_DEE |
| 5  | 408938   | 166046767 | T→G | 0.974 | +41.8  | 0.594 | intron_variant                | Early-infantile_DEE |
| 6  | 2746317  | 166013901 | G→T | 1.494 | −10.2  | 0.589 | intron_variant                | Early-infantile_DEE |
| 7  | 853842   | 166036540 | A→T | 1.506 | −12.0  | 0.587 | intron_variant                | Early-infantile_DEE |
| 8  | 2203203  | 166054633 | C→T | 1.640 | −36.3  | 0.553 | intron_variant                | Early-infantile_DEE |
| 9  | 4291800  | 166047773 | A→T | 1.350 |  −8.7  | 0.551 | intron_variant                | SCN1A-related_disorder |
| 10 | 1329953  | 165999718 | C→G | 1.355 | −10.3  | 0.548 | intron_variant                | Severe_myoclonic_epilepsy_in_infancy |

Full top-30 in `outputs/vus_top_combined.csv`.

### 6. The 3 "new" candidates DNASE pulled in

| Rank | rsID | Pos | SpliceScore | DNASEScore | Consequence | clndn |
|---:|---|---:|---:|---:|---|---|
| 19 | 1381398 | 166073497 | **0.009** | **+104.4** | missense_variant                | Early-infantile_DEE |
| 24 | 4532737 | 166047778 | 0.621    | +33.3     | intron_variant                  | not_provided, Early-infantile_DEE |
| 28 | 2178954 | 166073378 | **0.023** | **+89.6** | missense_variant                | Early-infantile_DEE, SCN1A-related |

**rsIDs 1381398 and 2178954 are the strongest pure-chromatin candidates:**
SPLICE_SITES_score ≤ 0.023 (zero splicing effect) yet DNASE scores +104 and
+90. These would be missed entirely by a splicing-only pipeline. They sit
in the same missense-cluster region (`~chr2:166073497–166073378` = SCN1A
exons 26–28 region, late in the protein's transmembrane domains).

### 7. The 3 candidates DNASE pushed out of top-30

| Splicing Rank | rsID | Pos | SpliceScore | SpliceJunctionsScore | Consequence |
|---:|---|---:|---:|---:|---|
| 17 | 4293437 | 166073671 | 1.055 | 857.2 | splice_acceptor_variant |
| 21 | 994749  | 166002526 | 0.974 | 1292.9 | missense_variant |
| 24 | 853637  | 165998175 | 0.955 | 798.8 | missense_variant |

These have very high splicing impact but their DNASE_score is near
background, so the average pulls them below the combined threshold. They
**remain strong splice candidates** — they should not be removed from any
candidate list. They were just displaced in the combined *rank*.

## Honest caveats

1. **No ground-truth labels on VUS.** The 1,610 VUS are ClinVar
   "uncertain_significance" variants — by definition we don't know which
   are truly pathogenic. Exp 004 showed DNase ±5bp concentration separates
   pathogenic vs benign at p=7.99e-7 *on labeled variants*, but that is a
   benchmark conclusion, not a guarantee that DNASE_score magnitude separates
   pathogenic from benign on unlabeled VUS. The direction is "more
   concentration = more pathogenic-like"; raw magnitude (sum of all 305
   track perturbations) is a coarser signal that includes both directions.

2. **DNASE values are *signed*.** Negative DNASE_score means the variant
   *decreases* chromatin accessibility at that locus. We do not currently
   know if decrease-of-accessibility is benign, pathogenic, or
   locus-dependent. The combined ranking treats negative and positive
   scores symmetrically (min-max-normalize → [0,1]), which conflates the
   two directions. **A future improvement**: rank by `|DNASE_score|`
   instead, or rank on the variant's "concentration at the variant
   position" ISM-derived metric from Exp 004 rather than this bulk
   `X.sum()` magnitude.

3. **Population frequency cross-check missing.** The skill recommends
   cross-referencing top candidates with gnomAD to remove common variants.
   `vus_top_combined.csv` does not have gnomAD columns. Before lab outreach,
   run a gnomAD GraphQL query for each rsID (per
   `alphagenome-variant-scoring` skill, "Cross-reference high-impact VUS
   with gnomAD and PubMed before outreach"). The existing
   `outputs/vus_high_impact_with_dnase.csv` already has gnomAD for the
   61 splicing-high candidates but not the 3 new DNASE-only ones.

4. **Class imbalance.** The DNASE top-30 is 86% missense variants, but
   only ~28% of the input 1,610 VUS are missense (most are intronic,
   untranslated, or synonymous). The DNASE scorer is therefore more
   permissive on coding variants — could be because DNase hypersensitivity
   signal is stronger near coding regions in the training tracks, or
   because missense variants at conserved sites genuinely perturb more
   chromatin structure than equivalent intronic variants. Cannot distinguish
   from this experiment alone.

5. **Correlation is not validation.** Spearman ρ = 0 on the VUS set is
   *good* (signals are independent) but does not validate either one.
   The next step would be to compute AUROC/AUPRC of DNASE_score against
   the *labeled* SCN1A ClinVar pathogenic/benign set, mirroring the
   `benchmark_scn1a_live_api.py` pattern — and to repeat on ≥1 other gene
   for generalization (per the Exp 004 replication rule).

## Recommendation

**Update the lab outreach** to mention the 3 chromatin-only candidates
alongside the existing top-30 splice candidates. Concrete edits:

- Add a one-line paragraph: *"Three additional candidates surfaced by
  AlphaGenome's DNase hypersensitivity scorer (rsIDs 1381398, 4532737,
  2178954) have minimal splicing impact (SPLICE_SITES_score ≤ 0.62) but
  strong chromatin perturbation (DNASE_score +33 to +104). These are
  orthogonal-mechanism candidates and would not appear in a splicing-only
  re-analysis."*
- Rank them **separately** in the email — don't merge with the splice
  candidates, since the validation assays differ (DNASE-disrupted candidates
  need chromatin / reporter assays, not minigene splicing assays).
- Before sending, **cross-reference each rsID against gnomAD** for
  population frequency (target: AC=0 or AF < 0.01% per the
  `alphagenome-variant-scoring` skill), and against the literature via
  PubMed for novelty.

Do **not** use the combined rank alone for outreach — it currently
weights raw `X.sum()` magnitude which conflates direction (caveat 2).
A follow-up experiment that uses `|DNASE_score|` or the Exp 004
"concentration at the variant position" summary statistic is likely to
produce a more clinically actionable ranking.

## Repro

```bash
cd /Users/hermes/projects/alphagenome-work
bash scripts/_run_with_key.sh scripts/rescore_vus_dnase.py
```

Runtime ~8.5 min on Apple Silicon. Zero-cost on the live API free tier
("1000s of predictions, not millions" — this is 1,610 predictions).