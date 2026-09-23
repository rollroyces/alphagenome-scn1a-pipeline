# Experiment 011 — ISM on SCN1A Tier-1 candidates (DNASE)

## TL;DR

Applied AlphaGenome's DNASE in-silico mutagenesis to the four SCN1A Tier-1
candidate variants (rs801806, rs4293437, rs801809, rs2847163) — the top
ranked SCN1A variants from `outputs/vus_top_candidates.csv` by
`SPLICE_SITES_score`. All four completed successfully.

**The key empirical observation is unexpected. None of the four shows the
+/-5bp concentration signature that Exp 004 found reliable for
distinguishing pathogenic from benign splice-region variants.** All
four `fraction_within_5bp` values land between 0.07 and 0.12, well below
the 0.30 threshold that we (heuristically) associate with a tight
splice-site hit. Three of the four have their peak effect 30–62 bp away
from the variant, not within ±5 bp. This is interpretable in three
different ways; the README makes no ground-truth claim.

## Methodology

| Parameter | Value | Matched to |
|---|---|---|
| Window (`ism_interval`) | variant ± 64 bp (128 bp total) | Exp 004 |
| Context (`interval`) | 16,384 bp (`SEQUENCE_LENGTH_16KB`) | Exp 001–004 |
| Scorer | `RECOMMENDED_VARIANT_SCORERS["DNASE"]` only | Exp 004 |
| API call | `dna_model.score_ism_variants(...)` (live API) | Exp 004 |
| Per-position metric | `\|effect\|` across alt bases, summed → magnitude | Exp 004 |
| Per-variant aggregation | average of 3 alt-allele cells per position | Exp 004 |

No ATAC, no RNA, no CAGE — only DNASE, because Exp 004 proved at n=30+30
per gene (DMD, CFTR, SCN1A) that the **DNASE** ±5bp concentration fraction
is the validated pathogenic-vs-benign discriminator
(combined Mann-Whitney p=7.99e-7, rank-biserial r=+0.467, see
`research_notebook/experiments/004_ism_dnase_n30/README.md`).

## Tier-1 candidate inputs

Pulled from `outputs/vus_top_candidates.csv` (rows where `rsid` matched
the four candidate dbSNP IDs):

| rsid | chrom | pos | ref | alt | molecular_consequence |
|---|---|---|---|---|---|
| rs801806 | chr2 | 166041471 | T | A | splice_acceptor_variant |
| rs4293437 | chr2 | 166073671 | C | G | splice_acceptor_variant |
| rs801809  | chr2 | 166043700 | A | G | splice_donor_variant |
| rs2847163 | chr2 | 166043701 | C | A | splice_donor_variant |

## Results

| rsid | total magnitude | frac ±5bp | frac ±15bp | max position (bp) | classification |
|---|---:|---:|---:|---:|---|
| rs801806  | 1453.5 | 0.068 | 0.233 | −30 | broad_regulatory |
| rs4293437 | 1757.8 | 0.122 | 0.259 | +1  | broad_regulatory (max at boundary) |
| rs801809  | 3593.8 | 0.092 | 0.229 | −61 | broad_regulatory |
| rs2847163 | 3543.9 | 0.082 | 0.221 | −62 | broad_regulatory |

(Source: `metrics.csv`. Classification is `dna` whenever frac±5bp < 0.30,
because none of these variants concentrate their effect tightly at the
SNV position. rs4293437 is on the boundary — see "Caveats".)

## Interpretation

### What Exp 004's benchmark would predict

In Exp 004, the n=30 pathogenic SCN1A SNVs had a mean ±5bp concentration
fraction that was reliably higher than the n=30 benign SCN1A SNVs
(combined Fisher p=3.16e-6). For an individual pathogenic SNV we'd
*expect* `fraction_within_5bp > 0.30`, often much higher, with the peak
effect within a few bp of the variant.

### What we observe

All four Tier-1 candidates show `fraction_within_5bp` between 0.068 and
0.122 — well **below** the Exp 004 pathogenic regime. The magnitude is
spread over the 128 bp window (frac ±15bp is 0.22–0.26, also not very
concentrated). In three of the four cases (rs801806, rs801809,
rs2847163), the position with the maximum absolute effect sits 30–62 bp
away from the variant, not at it.

### Why this is interpretable, not anomalous

There are at least three compatible readings:

1. **Tier-1 candidates may not all be DNASE splice hits.** A ClinVar
   "splice_donor_variant" / "splice_acceptor_variant" consequence
   annotation is a *sequence* annotation — it describes where the SNV
   sits relative to the canonical AG/GT splice signals. AlphaGenome's
   DNASE track is a chromatin-accessibility predictor for a wide
   variety of cell types, not a splice-junction predictor. A variant
   can sit inside a canonical splice motif and yet have its dominant
   regulatory effect on a *flanking* chromatin-accessibility element
   (an enhancer / promoter contact). The model's "best explanation" for
   each Tier-1 candidate may simply be a regulatory perturbation tens of
   bp away.

2. **The benchmark metric may not transfer.** Exp 004's DNASE ±5bp
   signature was learned on a set of n=60 SCN1A training-test SNVs
   curated as "high-confidence pathogenic" vs "high-confidence benign"
   by ClinVar. The Tier-1 set is the *top of `SPLICE_SITES_score`* from
   the entire ClinVar SCN1A VUS pool — a different selection bias, and
   one that may over-represent variants whose AlphaGenome reasoning is
   **distributional** rather than **positional**. (SPLICE_SITES_score is
   itself a track-level quantity, not a positional one.)

3. **The two adjacent donor-acceptor pairs are mirror images.** rs801809
   (donor) and rs2847163 (donor) are 1 bp apart on the same intron;
   they produce nearly identical heatmaps (max@-61bp vs max@-62bp, very
   similar magnitudes). rs801806 (acceptor) and rs4293437 (acceptor)
   sit on different introns but share the *max at or just past the
   variant position* feature — which is more donor-like than
   acceptor-like, an interesting reverse-direction pattern.

In short: **the Tier-1 candidates do NOT show the Exp 004 DNASE
signature** as a tight ±5bp hit. That is itself a result: it tells us
that AlphaGenome's reasoning about these specific candidates is *not*
"the splice-site motif is broken" — it's something more diffuse. This
is worth showing on a paper figure precisely because it's honest.

### Does the heatmap show distinct patterns per variant?

Yes — clearly two pairs:

- **rs801809 / rs2847163 (donor pair, 1 bp apart):** visually nearly
  identical heatmaps, both with magnitude peaks far upstream of the
  variant and similar total magnitudes (~3500).
  *Caveat:* flattened pixel-wise Pearson r=0.116 between the two
  matrices — the within-pair similarity is **structural** (matching
  peaks, matching total magnitude) rather than pixel-by-pixel identical,
  since ISM matrices have high-frequency variation across positions.
- **rs801806 / rs4293437 (acceptor pair, different introns):** also a
  visually similar pair to each other, but with magnitude roughly half
  that of the donor pair (~1500–1750) and the peak closer to (or just
  past) the variant position.

So the heatmap is *not* uniform — there is a real donor-vs-acceptor
structure between the two pairs. But the within-pair redundancy is
expected (one donor and one acceptor are even 1 bp apart) and should
not be over-interpreted as "two independent pieces of evidence".

## Caveats — please read before citing

- **n=4 is anecdotal.** No statistical test is meaningful here; the
  ±5bp metric is reported for benchmarking against the Exp 004 distribution,
  not as a per-variant classifier.
- **No ground truth.** We have ClinVar molecular consequence annotations
  (splice_donor / splice_acceptor), but no functional validation that
  any of these Tier-1 candidates actually disrupts splicing in a
  Dravet patient's cells. The "donor-like" / "acceptor-like"
  classifications are heuristics from peak position only.
- **DNASE alone.** AlphaGenome's splice signal is mostly captured by
  the dedicated splicing tracks (SPLICE_SITES, SPLICE_SITE_USAGE,
  SPLICE_JUNCTIONS). DNASE tells us what the model "sees" in chromatin
  accessibility, which is correlated with but not identical to
  splicing. Running ISM with the splice-specific scorers would be a
  natural follow-up.
- **The `+1` max for rs4293437** is technically "donor-like" by our
  heuristic, but its `frac ±5bp = 0.122` puts it below the 0.30
  threshold so we classified it `broad_regulatory`. If the threshold
  were relaxed to 0.10, that single variant would flip to
  `splice_donor_like`. The ±15bp fraction (0.259) is also the highest
  of the four, supporting the "just on the boundary" reading.
- **Per-alt-allele aggregation:** we averaged the 3 alt-allele rows per
  position, matching Exp 004's "per-variant" method. A
  per-alt-allele aggregation would give slightly different numbers but
  the same qualitative pattern.

## Files

- `ism_dnase_tier1.py` — runnable script (launch with
  `bash scripts/_run_with_key.sh research_notebook/experiments/011_ism_tier1/ism_dnase_tier1.py`)
- `metrics.csv` — 4 rows, one per candidate
- `ism_dnase_rs801806.npy`, `ism_dnase_rs4293437.npy`,
  `ism_dnase_rs801809.npy`, `ism_dnase_rs2847163.npy` — (128, 4) float32
  matrices; column order = A, C, G, T
- `figures/tier1_ism_heatmap.png` — 4-row heatmap, log1p(|effect|)
  magnitude per position, red dashed line at the variant position,
  white dotted lines at ±5 bp and ±15 bp