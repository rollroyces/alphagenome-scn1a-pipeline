# Experiment 006 — Tissue-specific vs averaged DNASE for SCN1A (brain) and DMD (muscle)

**Date:** September 23, 2026
**Status:** **MIXED / NEGATIVE** — Tissue-filtered DNASE does NOT consistently
outperform averaged DNASE. For SCN1A it gives a small win on max-abs; for DMD
it gives essentially no difference (averaged is marginally tighter on max-abs,
tissue-filtered is marginally better on mean-abs). The averaged 305-track
DNASE is already close to ceiling for these two genes at this n.

## Background and motivation

Experiment 004 (n=30+30 per gene, DNASE ISM ±5bp concentration) showed:
- DMD: p=7.97e-6, r=+0.695 (large effect, muscle gene)
- SCN1A: p=0.007, r=+0.421 (medium effect, brain gene)
- CFTR: p=0.318 (no effect, epithelial gene)

Experiment 005 (SPLICE_JUNCTIONS, brain-filtered vs averaged across 5 genes)
showed tissue-filtered tracks did NOT significantly outperform averaged tracks
(Wilcoxon on AUPRC deltas, p>0.05).

This experiment asks the analogous question for **DNASE**: does filtering the
305 DNASE tracks down to brain-relevant tracks (for SCN1A) or muscle-relevant
tracks (for DMD) yield a stronger pathogenic-vs-benign separation than using
all 305 tracks?

## Method

### Variant selection — IDENTICAL to Experiment 004

Top-30 pathogenic + bottom-30 benign SNVs per gene, ranked by
`SPLICE_SITES_score`, then per-variant aggregation (average across the 3 alt
alleles of each `(pos, ref)` group):

| Gene   | CSV                                            | n_path unique | n_ben unique |
|--------|------------------------------------------------|---------------|--------------|
| SCN1A  | `outputs/benchmark_scn1a_live_api_raw.csv`     | 20            | 28           |
| DMD    | `outputs/cross_disease_dmd_raw.csv`            | 24            | 29           |

CFTR was excluded (task brief did not include it; we already have its
Exp 004 null result, r=+0.088, p=0.318).

### Scoring — IMPORTANT deviation from the task brief

The task brief said: "Run DNASE CenterMaskScorer ISM with 64bp window each
side (same as Exp 004)". This is **not feasible** as written:

1. `dna_client.score_ism_variants` does NOT accept `ontology_terms` at the
   API level. Only `atlas.query_variants` supports tissue filtering.
2. The ISM AnnData in Exp 004 is pre-aggregated across all 305 DNASE tracks
   **before** the matrices are returned (the saved `.npy` files are
   `(128 positions, 4 bases)` — track-level resolution is permanently lost
   in the API response).

We therefore used **`atlas.query_variants`** (same endpoint as Exp 005) for
DNASE scoring. This returns AnnData of shape `(n_variants, 305)` with track
metadata in `ad.var` (biosample_name, ontology_curie, etc.). Per-variant
score = `max(|X[i, :]|)` across the track axis — exactly as the task brief
specified ("Per-variant score = max-abs across tracks (NOT junction ×
tracks — DNASE doesn't have junctions)").

This is the *only* way to do tissue filtering at the track level with the
current AlphaGenome API. Documented up-front so the comparison is
interpretable. The comparison still answers the substantive question:
"does tissue-filtered DNASE outperform averaged DNASE on pathogenic-vs-
benign separation?" — using identical variants and an identical
per-variant aggregation as Exp 004.

### Conditions

| Condition | `ontology_terms` | DNASE tracks |
|-----------|-------------------|---------------|
| `avg`     | (none)            | 305 (all)     |
| `tissue`  | brain ATS for SCN1A; muscle ATS for DMD | 20 (SCN1A), 15 (DMD) |

Tissue tracks were identified by regex on `biosample_name`:

- **Brain (SCN1A, 20 tracks):** biosample_name matches any of `brain | neuron |
  cortex | hippocamp | cerebell | neural | glia | spinal | forebrain |
  midbrain | dorsolateral` (case-insensitive).
- **Muscle (DMD, 15 tracks):** biosample_name matches any of `muscle | myocyte |
  skeletal | cardiac | myoblast`, AND NOT `brain vasculature` (which is a
  brain pericyte, not skeletal/cardiac muscle).

Full biosample lists for tissue tracks are printed by the script at run
time.

### Statistical test

Mann-Whitney U, `alternative='greater'` (pathogenic > benign predicted),
with effect size `r = 2U/(n1*n2) - 1` (Kerby 2014). Sign convention: `r>0`
means pathogenic stochastically greater than benign.

Reported metrics: **max-abs** (primary, matches brief) and **mean-abs**
(secondary, mirrors Exp 004's "concentration" metric at the track-axis
analogue).

## Results

### Per-variant max-abs (primary metric — task brief specified max-abs across tracks)

| Gene   | Cond   | n_path | n_ben | n_tracks | path mean | ben mean | U    | p-value | r      | Sig |
|--------|--------|--------|-------|----------|-----------|----------|------|---------|--------|-----|
| SCN1A  | avg    | 20     | 28    | 305      | 0.2162    | 0.1193   | 434  | 6.64e-4 | +0.550 | *** |
| SCN1A  | tissue | 20     | 28    | 20       | 0.1514    | 0.0669   | 458  | 1.03e-4 | +0.636 | *** |
| DMD    | avg    | 24     | 29    | 305      | 0.3448    | 0.0975   | 652  | 2.93e-8 | +0.874 | *** |
| DMD    | tissue | 24     | 29    | 15       | 0.2501    | 0.0558   | 651  | 3.24e-8 | +0.871 | *** |

### Per-variant mean-abs (secondary)

| Gene   | Cond   | path mean | ben mean | U    | p-value | r      | Sig |
|--------|--------|-----------|----------|------|---------|--------|-----|
| SCN1A  | avg    | 0.0503    | 0.0217   | 481  | 1.38e-5 | +0.718 | *** |
| SCN1A  | tissue | 0.0445    | 0.0208   | 479  | 1.65e-5 | +0.711 | *** |
| DMD    | avg    | 0.0600    | 0.0214   | 645  | 5.85e-8 | +0.853 | *** |
| DMD    | tissue | 0.1100    | 0.0216   | 658  | 1.60e-8 | +0.891 | *** |

### Comparison: tissue vs averaged

| Gene   | Metric    | Δp (tissue − avg, neg = tissue better) | Δr (tissue − avg, pos = tissue better) | Verdict |
|--------|-----------|---------------------------------------|----------------------------------------|---------|
| SCN1A  | max-abs   | **−5.6e-4** (tissue 6.5× smaller)     | **+0.086** (tissue larger)             | tissue **wins** |
| SCN1A  | mean-abs  | +2.8e-6 (tissue slightly worse)       | −0.007 (tissue slightly worse)         | essentially **tied** |
| DMD    | max-abs   | +3.1e-9 (tissue slightly worse)       | −0.003 (tissue slightly worse)         | essentially **tied** |
| DMD    | mean-abs  | **−4.3e-8** (tissue 3.7× smaller)     | **+0.038** (tissue larger)             | tissue **wins** |

## Headline interpretation

**Mixed / negative result.** The averaged (305-track) DNASE scoring already
separates pathogenic from benign very strongly on both genes (SCN1A
p≈7e-4, DMD p≈3e-8). Tissue filtering to ~15-20 tissue-specific tracks:

- **SCN1A (brain):** gives a modest improvement on max-abs (p drops
  6.5×, r gains +0.086). But mean-abs is essentially unchanged. The
  gain is small and not dramatic.
- **DMD (muscle):** gives essentially no improvement. The max-abs p
  is fractionally worse; mean-abs is fractionally better. Both are
  well within noise.

This mirrors the **Exp 005 SPLICE_JUNCTIONS** result (tissue-filtered
did not significantly outperform averaged across 5 genes). Together
these are **two independent scorers** (SPLICE_JUNCTIONS in Exp 005,
DNASE in Exp 006) reaching the same conclusion: **for pathogenic-vs-
benign classification in rare-disease Mendelian genes, tissue-specific
track filtering does NOT reliably outperform the full averaged track
panel.**

### Possible reasons

1. **Averaged tracks already include the tissue tracks.** Averaging
   across 305 tracks dilutes noise and retains signal — the relevant
   tissue tracks still drive the average up when they fire, while the
   280 non-tissue tracks act as a regularizer.
2. **Sample sizes are small (~20-29 per arm).** Effects need to be
   large to show up. DMD's averaged p≈3e-8 is already at floor.
3. **The pathogenic variants in the cross-disease CSVs are dominated
   by splice donor/acceptor variants**, which create large DNASE
   changes regardless of tissue. The signal may be tissue-agnostic
   for this specific variant selection.
4. **The AlphaGenome paper itself notes** that tissue-specific
   recapitulation remains a known limitation, so it would not be
   surprising if tissue-specific tracks do not yet outperform
   averaged tracks for variant interpretation.

### When tissue filtering *might* help (future work)

- Variants with **subtle, tissue-restricted** effects (e.g. regulatory
  variants affecting only one cell type). These wouldn't show up in
  our selection (top-N by SPLICE_SITES_score, which is tissue-agnostic).
- **Mosaic or low-effect-size** variants where the noise from
  280 non-tissue tracks swamps the signal.

## Output files

| File | Description |
|------|-------------|
| `summary_by_gene_condition.csv` | Per-gene × per-condition Mann-Whitney U, p, r |
| `per_variant_scn1a_avg.csv` | All 30 SCN1A pathogenic + 30 benign variants, averaged condition |
| `per_variant_scn1a_tissue.csv` | Same variants, tissue-filtered condition |
| `per_variant_dmd_avg.csv` | All DMD variants, averaged condition |
| `per_variant_dmd_tissue.csv` | Same variants, tissue-filtered condition |
| `tissue_vs_averaged.png` | Strip plot of max-abs scores per condition per gene |
| `run.log` | Full run output |

## Reproducibility

```bash
bash scripts/_run_with_key.sh \
    research_notebook/experiments/006_tissue_specific_dnase/tissue_dnase_experiment.py
```

Total runtime: ~25 seconds (4 API calls × 30 variants × 2 genes ×
2 conditions).

## Constraints honored

- ✅ Did not modify any existing files in prior experiments.
- ✅ Did not modify `paper/preprint.md`.
- ✅ Used the SAME 20+28 SCN1A and 24+29 DMD variants as Exp 004.
- ✅ Reused code patterns from Exp 004 (variant selection, per-variant
  3-alt aggregation, Mann-Whitney + Kerby r) and Exp 005 (atlas client
  + ontology_terms filtering + biosample regex).
- ✅ Honest negative result documented (this README).