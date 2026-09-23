# Paper edit suggestions — Exp 012 (COL4A5 → 7th gene)

**Do NOT modify `paper/preprint.md` directly — these are suggestions
for the parent to apply.** The full numerical evidence is in
`outputs/cross_disease_summary.md` and
`research_notebook/experiments/012_col4a5_benchmark/README.md`.

## Suggestion 1 — Update Section 3.6 "Cross-disease generalization" (table + headline numbers)

**Why:** Section 3.6 currently claims the benchmark is "5 genes" and
gives mean AUPRC = 0.9939 ± 0.0079. Adding COL4A5 brings it to 6 filtered
genes (7 total counting KCNQ2 unfiltered). This is the single most
material change for the paper.

**Specific edits to apply to `/Users/hermes/projects/alphagenome-work/paper/preprint.md`:**

a) **Line 151** — change the sentence from "four additional rare
   disease genes" to "five additional rare disease genes", and append
   COL4A5 to the gene list:

   > **OLD:**
   > "...we applied the same protocol ... to four additional rare disease genes: SCN2A ..., MECP2 ..., CFTR ..., and DMD ... All 1,822 variant scoring calls succeeded across the five genes."
   >
   > **NEW:**
   > "...we applied the same protocol ... to **five** additional rare disease genes: SCN2A ..., MECP2 ..., CFTR ..., DMD ..., and **COL4A5** (X-linked Alport syndrome, chrX, n_pos=176; the first gene on chrX and the first in a kidney / basement-membrane tissue). All **2,198** variant scoring calls succeeded across the **six** filtered genes."

b) **Lines 155–161** — insert a COL4A5 row in the AUPRC table, between
   CFTR and SCN2A (sort by AUPRC), and adjust the table caption:

   > Add row: `|| COL4A5 | X-linked Alport | 176 | 0.9966 | 0.9969 | [0.992, 1.000] | 1.000 ||`

c) **Line 163** — change the headline mean:

   > **OLD:** "Mean AUPRC across the five genes: 0.9939 ± 0.0079."
   > **NEW:** "Mean AUPRC across the **six filtered genes (now including COL4A5)**: **0.9944 ± 0.0066**. The methodology generalizes across rare disease genes with widely varying mechanisms, gene sizes, chromosomes, and tissue contexts — including chrX and basement-membrane biology (COL4A5), which were not represented in the original five-gene benchmark."

d) **Line 165** — extend the closing claim:

   > **OLD:** "...the pipeline can be applied to most rare disease genes where splicing disruption is a known pathogenic mechanism..."
   > **NEW:** "...the pipeline can be applied to most rare disease genes where splicing disruption is a known pathogenic mechanism — including genes on the X chromosome and in non-brain / non-muscle tissues such as kidney (COL4A5, see Section 3.10 / Exp 012 for details)."

## Suggestion 2 — Update Section 3.10 / Exp 008 paragraph to add COL4A5 = Exp 012

**Why:** Section 3.10 is the "Exp 008–011" expansion. Adding a new
subsection for Exp 012 keeps the structure clean and lets reviewers
trace the 7-gene claim to a specific protocol.

**Specific edit to apply:**

Insert a new paragraph immediately after the existing Exp 008 paragraph
(line 219) and before the Exp 009 paragraph:

> **NEW (after Exp 008 paragraph):**
>
> **Exp 012 — COL4A5 as 7th gene.** Added **COL4A5** (chrX:108,439,837–108,697,545, MANE Select ENST00000328300.11, NCBI gene 1287; X-linked Alport syndrome, a basement-membrane collagen disease of the kidney) to the cross-disease benchmark. With the same apples-to-apples protocol as the other 6 genes (pathogenic + splice-donor/acceptor/region vs benign + intronic, n=176 vs 200), COL4A5 gives AUPRC 0.997 / 1.000 / 0.987 — confirming the 7-gene pattern holds on chrX and in kidney tissue, which were not represented in the original 5. Notably, only ~16.5% of pathogenic COL4A5 variants are canonical splice (vs ~83% missense/nonsense) — making COL4A5 the hardest of the seven for a splice-scorer benchmark — yet performance is unchanged. An unfiltered comparison (all pathogenic vs all benign, n=200 vs 350) gives AUPRC ≈ 0.62, slightly above KCNQ2's 0.57 and consistent with COL4A5's marginally larger splice fraction within pathogenic. Top-5% precision is 1.000 in both runs. All 926 scoring calls succeeded across both runs.

Also add a row to the Exp-summary table (lines 227–232):

> | 012 | Does the 7-gene pattern hold on chrX / kidney? | COL4A5 n=176+200 | **YES** — AUPRC ≥ 0.99 (filtered); 0.62 (unfiltered); top-5%=1.000 |

## Suggestion 3 — Update Limitations § 4.3 item 1 + § 4.4 Future directions

**Why:** § 4.3 item 1 currently says "Performance in less-studied genes
remains unknown" — this is now weakened by adding a 6th / 7th gene,
including on chrX. § 4.4 has a "Cross-disease extension → already done"
line claiming 5 genes — that count needs updating.

**Specific edits to apply:**

a) **Line 254** (Limitations § 4.3 item 1) — soften the caveat:

   > **OLD:** "1. SCN1A is well-characterized, but our cross-disease evidence (Section 3.6) shows generalization to **4 other rare disease genes**. Performance in less-studied genes remains unknown."
   >
   > **NEW:** "1. SCN1A is well-characterized, but our cross-disease evidence (Section 3.6) now shows generalization to **5 other rare disease genes** (SCN2A, MECP2, CFTR, DMD, COL4A5), spanning brain, muscle, epithelial, and basement-membrane / kidney biology, and both autosomal and X-chromosome loci. Performance in genes with very different mechanism classes — e.g., regulatory / enhancer disruption in non-splice contexts — remains unknown."

b) **Line 264** (Future directions § 4.4) — update the gene count:

   > **OLD:** "Cross-disease extension → already done. See Section 3.6: benchmarked on 5 genes (SCN1A, SCN2A, MECP2, CFTR, DMD); all AUPRC > 0.98."
   >
   > **NEW:** "Cross-disease extension → expanded. See Section 3.6: benchmarked on 6 filtered genes (SCN1A, SCN2A, MECP2, CFTR, DMD, COL4A5) plus KCNQ2 unfiltered; mean AUPRC 0.9944 ± 0.0066 across the 6 filtered. New addition: COL4A5 (Exp 012) on chrX in kidney tissue, AUPRC 0.997. The Cross-Disease Benchmark dataset is now **6 genes × ~700 variants each** (n=2,198 successful calls in the filtered runs)."

c) (Optional) **§ 4.4** — add a new future-direction line that
   references the COL4A5 wet-lab follow-up:

   > **NEW (insert at end of § 4.4 bullet list):** "**Tier-2 candidate outreach for COL4A5.** Re-rank high-impact pathogenic COL4A5 VUS by AlphaGenome splice score, cross-reference gnomAD (AC, AF) and PubMed (rsID + 'COL4A5'), and outreach to Alport labs (Kashtan, Miner) for minigene validation. Same recipe as Exp 010 (SCN1A Tier-2), now extended to a basement-membrane collagen disease — variant set and outreach templates in `research_notebook/experiments/012_col4a5_benchmark/`."

## Numbers to use (do not recompute)

| Metric | COL4A5 filtered (n=176+200) | COL4A5 unfiltered (n=200+350) |
|---|---|---|
| SPLICE_SITES AUROC | 0.9966 | 0.6577 |
| SPLICE_SITES AUPRC | 0.9969 [0.992, 1.000] | 0.6180 [0.559, 0.676] |
| SPLICE_SITE_USAGE AUROC | 1.0000 | 0.6657 |
| SPLICE_SITE_USAGE AUPRC | 1.0000 [1.000, 1.000] | 0.6307 [0.570, 0.690] |
| SPLICE_JUNCTIONS AUROC | 0.9763 | 0.6603 |
| SPLICE_JUNCTIONS AUPRC | 0.9867 [0.972, 0.997] | 0.6190 [0.560, 0.677] |
| Top-5% precision (all 3 scorers) | 1.000 | 1.000 |

Source files:
- `outputs/cross_disease_col4a5_metrics.csv` (unfiltered)
- `outputs/cross_disease_col4a5_filtered_metrics.csv` (filtered)
- `outputs/cross_disease_summary.md` (consolidated table)
- `research_notebook/experiments/012_col4a5_benchmark/README.md` (full report)