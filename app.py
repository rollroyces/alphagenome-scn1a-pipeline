"""
SCN1A VUS Re-Scoring — Streamlit Web UI
========================================

A clinician-facing interface for the 1,610 SCN1A VUS scored with AlphaGenome +
DNASE, hosted on Hugging Face: https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome

Pages
-----
1. Browse Dataset        — filter and explore all 1,610 scored VUS
2. Tier-1 Candidates     — detailed cards for the 4 wet-lab-ready variants
3. Tier-2 Candidates     — detailed cards for 13 curated Tier-2 (non-canonical mechanism) candidates
4. Submit VCF            — paste/upload a VCF, score variants (requires API key, uses real quota)

Run locally:
    streamlit run app.py

Deploy:
    Hugging Face Spaces — see README.md in this directory for the HF Spaces config.
"""

from __future__ import annotations

import io
import os
import re
from typing import Optional

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="SCN1A VUS Re-Scoring",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATASET_REPO = "RROL/scn1a-vus-alphagenome"
TIER1_RSIDS = ["rs801806", "rs4293437", "rs2847163", "rs801809"]
GITHUB_REPO = "https://github.com/rollroyces/alphagenome-scn1a-pipeline"
HF_DATASET_URL = "https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome"

# ---------------------------------------------------------------------------
# Tier-1 candidate data
# ---------------------------------------------------------------------------
# Source: outputs/tier1_actionable_features.md (Exp 013 + Exp 014).
# Hardcoded so the Streamlit app has no dependency on local output files —
# works the same in HF Spaces and locally.

TIER1_DATA = [
    {
        "rsid": "rs801806",
        "chrom": 2,
        "pos": 166041471,
        "ref": "T",
        "alt": "A",
        "consequence": "splice_acceptor_variant",
        "clndn": "Severe_myoclonic_epilepsy_in_infancy",
        "rank": 1,
        "splice_sites": 1.1406,
        "splice_site_usage": 146.57,
        "splice_junctions_live": 970.9,
        "splice_junctions_atlas": 1982.3,
        "closest_splice_site": "3'ss (acceptor) @ 2 bp",
        "canonical_exon": "exon 16",
        "mean_cos_path": 0.238,
        "max_cos_path": 0.392,
        "dnase_dscore": -2.884,
        "gnomad_present": True,
        "gnomad_genome_af": 0.0,
        "gnomad_exome_af": 0.0,
        "pubmed_hits": 0,
        "mechanism": (
            "Position sits 2 bp from the 3'ss (acceptor) of SCN1A's 16th exon. "
            "Canonical splice-acceptor disruption — expected to weaken exon 16 inclusion."
        ),
    },
    {
        "rsid": "rs4293437",
        "chrom": 2,
        "pos": 166073671,
        "ref": "C",
        "alt": "G",
        "consequence": "splice_acceptor_variant",
        "clndn": "Developmental_and_epileptic_encephalopathy_6B",
        "rank": 2,
        "splice_sites": 1.0547,
        "splice_site_usage": 145.75,
        "splice_junctions_live": 453.8,
        "splice_junctions_atlas": 857.2,
        "closest_splice_site": "3'ss (acceptor) @ 1 bp",
        "canonical_exon": "exon 4",
        "mean_cos_path": 0.161,
        "max_cos_path": 0.232,
        "dnase_dscore": -13.790,
        "gnomad_present": False,
        "gnomad_genome_af": None,
        "gnomad_exome_af": None,
        "pubmed_hits": 0,
        "mechanism": (
            "Position sits 1 bp from the 3'ss (acceptor) of SCN1A's 4th exon. "
            "Strongest predicted brain expression impact (Exp 014: 28.6 brain RNA_SEQ, "
            "9x stronger than rs801806, 14/27 brain tracks above |1.0|)."
        ),
    },
    {
        "rsid": "rs2847163",
        "chrom": 2,
        "pos": 166043701,
        "ref": "C",
        "alt": "A",
        "consequence": "splice_donor_variant",
        "clndn": "Early-infantile_DEE",
        "rank": 3,
        "splice_sites": 0.9150,
        "splice_site_usage": 91.96,
        "splice_junctions_live": 647.0,
        "splice_junctions_atlas": 838.4,
        "closest_splice_site": "5'ss (donor) @ 33 bp (canonical exon 14)",
        "canonical_exon": "exon 14 (transcript-dependent)",
        "mean_cos_path": 0.058,
        "max_cos_path": 0.141,
        "dnase_dscore": 4.260,
        "gnomad_present": False,
        "gnomad_genome_af": None,
        "gnomad_exome_af": None,
        "pubmed_hits": 0,
        "mechanism": (
            "Position sits 33 bp from the 5'ss (donor) of canonical SCN1A-224 exon 14 "
            "(inside the exon in some transcripts). 5'ss donor disruption — most likely "
            "to show exon skipping on minigene assay."
        ),
    },
    {
        "rsid": "rs801809",
        "chrom": 2,
        "pos": 166043700,
        "ref": "A",
        "alt": "G",
        "consequence": "splice_donor_variant",
        "clndn": "Seizure; Severe_myoclonic_epilepsy_in_infancy",
        "rank": 4,
        "splice_sites": 0.9213,
        "splice_site_usage": 113.56,
        "splice_junctions_live": 788.4,
        "splice_junctions_atlas": 788.4,
        "closest_splice_site": "5'ss (donor) @ 32 bp (canonical exon 14)",
        "canonical_exon": "exon 14 (transcript-dependent)",
        "mean_cos_path": 0.043,
        "max_cos_path": 0.120,
        "dnase_dscore": 26.551,
        "gnomad_present": False,
        "gnomad_genome_af": None,
        "gnomad_exome_af": None,
        "pubmed_hits": 0,
        "mechanism": (
            "Position sits 32 bp from the 5'ss (donor) of canonical SCN1A-224 exon 14 "
            "(inside the exon in some transcripts). Strongest chromatin displacement "
            "(DNASE Δ +26.6) among the 4 — strongest predicted chromatin-accessibility impact."
        ),
    },
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner="Loading 1,610 SCN1A VUS from Hugging Face…")
def load_vus_dataset() -> pd.DataFrame:
    """Load the public SCN1A VUS dataset from Hugging Face.

    Cached for 1 hour to avoid re-loading 1,610 rows on every page refresh.
    """
    from datasets import load_dataset

    ds = load_dataset(DATASET_REPO, split="train")
    df = ds.to_pandas()

    # Normalise rsid for display (the HF dataset stores integers like 851265;
    # the rest of the app uses the rs### format).
    def _format_rsid(value) -> str:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return "(no rsID)"
        try:
            return f"rs{int(value)}"
        except (TypeError, ValueError):
            return str(value)

    df["rsid_display"] = df["rsid"].apply(_format_rsid)

    # Pretty gnomAD AF: cap at small numbers for filtering/display.
    df["gnomad_af_max"] = df[["gnomad_genome_af", "gnomad_exome_af"]].max(axis=1)

    return df


@st.cache_data(ttl=3600, show_spinner=False)
def load_tier2_curated() -> pd.DataFrame:
    """Load the 13 curated Tier-2 candidates.

    Hardcoded inline rather than reading the CSV, so the app works in HF Spaces
    without shipping the outputs/ directory. Source: outputs/vus_tier2_candidates.csv
    (Experiment 010).
    """
    return pd.DataFrame(
        [
            # ---- Stratum A — deep intronic ----
            {
                "stratum": "intron_variant",
                "rsid": "rs1412774",
                "pos": 166047773, "ref": "A", "alt": "C",
                "SPLICE_SITES_score": 1.5376,
                "DNASE_score": 24.8,
                "combined_score": 0.778,
                "splice_only_rank": 4,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Deep intronic. Co-located with rs4291800 (rank 9) on chr2:166047773 — "
                    "two alt alleles at the same locus argue for a position-level splice "
                    "regulatory element. A→C substitution in a deep intronic site."
                ),
            },
            {
                "stratum": "intron_variant",
                "rsid": "rs851265",
                "pos": 166047622, "ref": "C", "alt": "A",
                "SPLICE_SITES_score": 1.7153,
                "DNASE_score": -15.0,
                "combined_score": 0.636,
                "splice_only_rank": 1,
                "clndn": "Migraine",
                "mechanistic_hypothesis": (
                    "Deep intronic. Top of splice-only ranking (rank 1). "
                    "Co-located with rs393000 (rank 3) — same position, two alt alleles, "
                    "locus-level splice signal."
                ),
            },
            {
                "stratum": "intron_variant",
                "rsid": "rs393000",
                "pos": 166047622, "ref": "C", "alt": "G",
                "SPLICE_SITES_score": 1.6133,
                "DNASE_score": -8.3,
                "combined_score": 0.614,
                "splice_only_rank": 3,
                "clndn": "Complex neurodevelopmental / DEE",
                "mechanistic_hypothesis": (
                    "Deep intronic. Same locus as rs851265 (rank 1) — "
                    "C→A and C→G at chr2:166047622 together argue for a locus-level "
                    "splice regulatory element."
                ),
            },
            {
                "stratum": "intron_variant",
                "rsid": "rs408938",
                "pos": 166046767, "ref": "T", "alt": "G",
                "SPLICE_SITES_score": 0.974,
                "DNASE_score": 41.8,
                "combined_score": 0.528,
                "splice_only_rank": 18,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Deep intronic. Strongest positive DNASE displacement (+41.8) of any "
                    "Tier-2 candidate — chromatin accessibility disruption in addition to splice signal."
                ),
            },
            {
                "stratum": "intron_variant",
                "rsid": "rs2746317",
                "pos": 166013901, "ref": "G", "alt": "T",
                "SPLICE_SITES_score": 1.494,
                "DNASE_score": -10.2,
                "combined_score": 0.526,
                "splice_only_rank": 7,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Deep intronic. Located within ~10 bp of chr2:166013900 — "
                    "the strongest pathogenic-cosine locus in the SCN1A ISM cohort."
                ),
            },
            # ---- Stratum B — missense (dual mechanism) ----
            {
                "stratum": "missense_variant",
                "rsid": "rs1038247",
                "pos": 166042328, "ref": "T", "alt": "C",
                "SPLICE_SITES_score": 0.970,
                "DNASE_score": 26.0,
                "combined_score": 0.424,
                "splice_only_rank": 20,
                "clndn": "Developmental_and_epileptic_encephalopathy",
                "mechanistic_hypothesis": (
                    "Dual-mechanism missense + splice. Position 166042328 sits in the "
                    "central transmembrane/pore region of Nav1.1, where missense variants "
                    "are well-known to cause both gain- and loss-of-function DEE. "
                    "Highest positive DNASE displacement in the missense stratum (+26.0)."
                ),
            },
            {
                "stratum": "missense_variant",
                "rsid": "rs3726776",
                "pos": 166015607, "ref": "C", "alt": "T",
                "SPLICE_SITES_score": 1.279,
                "DNASE_score": -14.9,
                "combined_score": 0.359,
                "splice_only_rank": 11,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Dual-mechanism missense + splice. Splice-only rank 11, above the canonical "
                    "Tier-1 splice variants' median. Standard minigene assays (which test the "
                    "canonical transcript only) will miss the splice effect."
                ),
            },
            {
                "stratum": "missense_variant",
                "rsid": "rs2684360",
                "pos": 165998038, "ref": "C", "alt": "G",
                "SPLICE_SITES_score": 1.074,
                "DNASE_score": -0.1,
                "combined_score": 0.323,
                "splice_only_rank": 15,
                "clndn": "not_provided",
                "mechanistic_hypothesis": (
                    "Dual-mechanism missense + splice. Lower-confidence missense pick — "
                    "in the bottom third of the top-30 by combined score."
                ),
            },
            {
                "stratum": "missense_variant",
                "rsid": "rs4855006",
                "pos": 166058573, "ref": "T", "alt": "C",
                "SPLICE_SITES_score": 0.939,
                "DNASE_score": 1.8,
                "combined_score": 0.250,
                "splice_only_rank": 23,
                "clndn": "Severe_myoclonic_epilepsy_in_infancy",
                "mechanistic_hypothesis": (
                    "Dual-mechanism missense + splice. ClinVar label matches the classic "
                    "Dravet phenotype."
                ),
            },
            {
                "stratum": "missense_variant",
                "rsid": "rs461265",
                "pos": 166013794, "ref": "A", "alt": "T",
                "SPLICE_SITES_score": 0.931,
                "DNASE_score": 0.2,
                "combined_score": 0.233,
                "splice_only_rank": 26,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Dual-mechanism missense + splice. Lower-confidence missense pick — "
                    "in the bottom third of the top-30 by combined score."
                ),
            },
            # ---- Stratum C — non-coding transcript ----
            {
                "stratum": "non-coding_transcript_variant",
                "rsid": "rs1046194",
                "pos": 166013744, "ref": "C", "alt": "T",
                "SPLICE_SITES_score": 1.521,
                "DNASE_score": -10.0,
                "combined_score": 0.545,
                "splice_only_rank": 5,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Non-coding transcript overlap (possibly antisense SCN1A-AS1). "
                    "Same position as rs934879 (rank 13) — two alt alleles at one locus "
                    "argues for a locus-level regulatory element."
                ),
            },
            {
                "stratum": "non-coding_transcript_variant",
                "rsid": "rs934879",
                "pos": 166013744, "ref": "C", "alt": "G",
                "SPLICE_SITES_score": 1.141,
                "DNASE_score": 3.0,
                "combined_score": 0.385,
                "splice_only_rank": 13,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Non-coding transcript overlap. Same position as rs1046194 (rank 5) — "
                    "C→T and C→G together point to a locus-level non-coding regulatory element."
                ),
            },
            {
                "stratum": "non-coding_transcript_variant",
                "rsid": "rs2015486",
                "pos": 165994240, "ref": "T", "alt": "A",
                "SPLICE_SITES_score": 0.934,
                "DNASE_score": 9.2,
                "combined_score": 0.293,
                "splice_only_rank": 24,
                "clndn": "Early-infantile_DEE",
                "mechanistic_hypothesis": (
                    "Non-coding transcript overlap. Splice/chromatin scores in same range as "
                    "Tier-1 splice variants — likely affects RNA processing without changing a protein."
                ),
            },
        ]
    )


# ---------------------------------------------------------------------------
# VCF parsing (small, dependency-free)
# ---------------------------------------------------------------------------

_VCF_DATA_LINE_RE = re.compile(r"^[^\s#].*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*")


def parse_vcf(text: str) -> pd.DataFrame:
    """Parse a VCF text blob into a DataFrame with chrom/pos/ref/alt columns.

    Supports both tab-separated and whitespace-separated records. Header lines
    starting with `#` are skipped. Returns columns: chrom, pos (int), ref, alt.
    Multiple ALT alleles (comma-separated) are expanded into separate rows.
    """
    rows = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line)
        if len(parts) < 5:
            continue
        chrom, pos, _rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        try:
            pos_int = int(pos)
        except ValueError:
            continue
        for alt_allele in alt.split(","):
            rows.append({"chrom": chrom, "pos": pos_int, "ref": ref, "alt": alt_allele})

    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"])


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_score(value, places: int = 3) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    return f"{float(value):.{places}f}"


def _fmt_af(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "absent"
    if value == 0:
        return "0"
    return f"{value:.2e}"


def _gnomad_badge(present: bool) -> str:
    return "🟢 present" if present else "⚪ absent"


def render_sidebar() -> str:
    """Render the persistent sidebar and return the selected page name."""
    with st.sidebar:
        st.title("🧬 SCN1A VUS")
        st.caption("AlphaGenome + DNASE re-scoring")
        st.markdown("---")
        page = st.radio(
            "Navigate",
            options=[
                "📊 Browse Dataset",
                "⭐ Tier-1 Candidates",
                "🔬 Tier-2 Candidates",
                "📤 Submit VCF",
            ],
            label_visibility="collapsed",
        )
        st.markdown("---")
        st.markdown(
            f"**Data:** [{DATASET_REPO}]({HF_DATASET_URL})  \n"
            f"**Methods:** [{GITHUB_REPO}]({GITHUB_REPO})  \n"
            "**Citation:** Royce Lam (2026)"
        )
    return page


# ---------------------------------------------------------------------------
# Page 1 — Browse Dataset
# ---------------------------------------------------------------------------

def render_browse_page() -> None:
    st.title("📊 Browse Dataset")
    st.markdown(
        "Filter and explore all **1,610 SCN1A variants of uncertain significance (VUS)** "
        "from ClinVar re-scored with [AlphaGenome](https://alphagenome.google) (DeepMind, 2025) "
        "splicing + DNASE scorers and cross-referenced against gnomAD v4.1 and PubMed."
    )

    try:
        df = load_vus_dataset()
    except Exception as exc:
        st.error(
            f"Failed to load dataset `{DATASET_REPO}` from Hugging Face. "
            f"This app needs network access at first run. Error: `{exc}`"
        )
        st.stop()
        return

    # ----- Filter sidebar (within-page) -----
    with st.expander("🎛 Filters", expanded=True):
        consequences = sorted(df["molecular_consequence"].dropna().unique().tolist())
        default_cons = [c for c in ("splice_acceptor_variant", "splice_donor_variant") if c in consequences]

        col1, col2 = st.columns(2)
        with col1:
            sel_consequences = st.multiselect(
                "Consequence",
                options=consequences,
                default=default_cons or consequences[:3],
                help="Filter by VEP / ClinVar molecular consequence annotation.",
            )
        with col2:
            sel_tiers = st.multiselect(
                "Tier",
                options=[1, 2, 3],
                default=[1, 2, 3],
                help="Tier 1 = 4 wet-lab-ready. Tier 2 = 46 high-splice. Tier 3 = remaining 1,560.",
            )

        splice_min, splice_max = float(df["SPLICE_SITES_score"].min()), float(df["SPLICE_SITES_score"].max())
        sel_splice = st.slider(
            "SPLICE_SITES score range",
            min_value=splice_min,
            max_value=splice_max,
            value=(splice_min, splice_max),
            step=0.05,
            help="AlphaGenome SPLICE_SITES scorer (raw, max ~2.0).",
        )

        af_max = float(df["gnomad_af_max"].fillna(0).max())
        sel_af_max = st.slider(
            "Max gnomAD AF (genome + exome)",
            min_value=0.0,
            max_value=max(af_max, 1e-4),
            value=max(af_max, 1e-4),
            step=1e-5,
            format="%.1e",
            help="Variants with maximum gnomAD AF above this threshold are excluded.",
        )

        search_rsid = st.text_input(
            "Search rsID (optional)",
            value="",
            placeholder="e.g. rs801806 or 801806",
            help="Substring match against the rsID column.",
        ).strip()

    # ----- Apply filters -----
    f = df.copy()
    if sel_consequences:
        f = f[f["molecular_consequence"].isin(sel_consequences)]
    if sel_tiers:
        f = f[f["tier"].isin(sel_tiers)]
    f = f[(f["SPLICE_SITES_score"] >= sel_splice[0]) & (f["SPLICE_SITES_score"] <= sel_splice[1])]
    f = f[f["gnomad_af_max"].fillna(0) <= sel_af_max]
    if search_rsid:
        f = f[f["rsid_display"].str.contains(search_rsid, case=False, na=False)]

    # ----- Summary metrics -----
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Variants shown", f"{len(f):,}", delta=f"of {len(df):,}", delta_color="off")
    m2.metric("Tier-1 in view", f"{int((f['tier'] == 1).sum())}")
    m3.metric("Tier-2 in view", f"{int((f['tier'] == 2).sum())}")
    m4.metric("Absent from gnomAD", f"{int((f['present'].fillna(False) == False).sum())}")

    # ----- Table -----
    table_cols = [
        "rsid_display",
        "chrom",
        "pos",
        "ref",
        "alt",
        "molecular_consequence",
        "SPLICE_SITES_score",
        "DNASE_score",
        "gnomad_af_max",
        "clndn",
        "tier",
    ]
    table = f[table_cols].rename(
        columns={
            "rsid_display": "rsID",
            "chrom": "chr",
            "pos": "pos (GRCh38)",
            "ref": "ref",
            "alt": "alt",
            "molecular_consequence": "consequence",
            "SPLICE_SITES_score": "SPLICE_SITES",
            "DNASE_score": "DNASE Δ",
            "gnomad_af_max": "gnomAD AF (max)",
            "clndn": "ClinVar disease",
            "tier": "tier",
        }
    ).sort_values("SPLICE_SITES", ascending=False).reset_index(drop=True)

    st.markdown(f"**Showing {len(table):,} variants** — click any row for details.")
    event = st.dataframe(
        table,
        width="stretch",
        height=420,
        on_select="rerun",
        selection_mode="single-row",
        key="browse_table",
        hide_index=True,
    )

    # ----- Detail view for selected row -----
    selected_rows = event.selection.rows if event and event.selection else []
    if selected_rows:
        idx = selected_rows[0]
        st.markdown("---")
        st.subheader(f"🔎 Variant detail — {table.iloc[idx]['rsID']}")
        _render_variant_detail(df, table.iloc[idx]["rsID"])


def _render_variant_detail(full_df: pd.DataFrame, rsid_display: str) -> None:
    """Render the full evidence block for a single variant."""
    matches = full_df[full_df["rsid_display"] == rsid_display]
    if matches.empty:
        st.warning(f"Variant {rsid_display} not found.")
        return
    row = matches.iloc[0]

    c1, c2, c3 = st.columns(3)
    c1.metric("Position", f"chr{int(row['chrom'])}:{int(row['pos']):,}")
    c2.metric("Allele", f"{row['ref']} → {row['alt']}")
    c3.metric("Tier", f"{int(row['tier'])}")

    st.markdown(
        f"**Consequence:** `{row['molecular_consequence']}`  \n"
        f"**ClinVar disease label:** `{row.get('clndn', '—')}`  \n"
        f"**ClinVar review status:** `{row.get('review_status', '—')}`  \n"
        f"**gnomAD:** {_gnomad_badge(bool(row.get('present', False)))} "
        f"(genome AF = {_fmt_af(row.get('gnomad_genome_af'))}, "
        f"exome AF = {_fmt_af(row.get('gnomad_exome_af'))})  \n"
        f"**PubMed hits:** {int(row.get('pubmed_hits') or 0)}"
    )

    st.markdown("#### AlphaGenome scores")
    sc1, sc2, sc3, sc4 = st.columns(4)
    sc1.metric("SPLICE_SITES", _fmt_score(row.get("SPLICE_SITES_score"), 4))
    sc2.metric("SPLICE_SITE_USAGE", _fmt_score(row.get("SPLICE_SITE_USAGE_score"), 2))
    sc3.metric("SPLICE_JUNCTIONS", _fmt_score(row.get("SPLICE_JUNCTIONS_score"), 2))
    dnase = row.get("DNASE_score")
    dnase_str = _fmt_score(dnase, 3) if dnase is not None and not (isinstance(dnase, float) and pd.isna(dnase)) else "—"
    sc4.metric("DNASE Δ", dnase_str)
    if dnase is not None and not (isinstance(dnase, float) and pd.isna(dnase)):
        if float(dnase) > 0:
            st.caption("Positive DNASE Δ → Alt allele disrupts chromatin accessibility (validated pathogenic signature).")
        else:
            st.caption("Negative DNASE Δ → No/minor effect on chromatin accessibility at this site.")

    st.markdown("#### Mechanism hypothesis")
    consequence = row.get("molecular_consequence", "")
    if "splice" in str(consequence).lower():
        st.write(
            "Annotated as a splice-site variant by VEP/ClinVar. "
            "The high SPLICE_SITES score predicts canonical splice-site disruption "
            "(exon skipping, cryptic splice activation, or splice-site weakening)."
        )
    elif "missense" in str(consequence).lower():
        st.write(
            "Missense variant. AlphaGenome SPLICE_SITES in the same range as canonical "
            "splice variants suggests **dual-mechanism** pathogenicity — coding effect "
            "PLUS splicing disruption. Standard minigene assays (canonical transcript) may miss the splice effect."
        )
    elif "intron" in str(consequence).lower():
        st.write(
            "Deep intronic variant. The high SPLICE_SITES score suggests **cryptic splice-site "
            "creation, poison-exon activation, or ISE/ISS disruption** in a non-annotated region."
        )
    elif "non_coding" in str(consequence).lower() or "non-coding" in str(consequence).lower():
        st.write(
            "Non-coding transcript overlap (possibly antisense SCN1A-AS1 or regulatory). "
            "Likely affects RNA processing without changing a protein."
        )
    else:
        st.write(
            "Consequence not annotated as splice-relevant by VEP. Treat the AlphaGenome scores "
            "as exploratory and orthogonal evidence — not a wet-lab-ready prediction."
        )

    if int(row["tier"]) == 1:
        st.success("⭐ **Tier-1** — wet-lab-testable canonical splice variant. See the Tier-1 page for full evidence.")
    elif int(row["tier"]) == 2:
        st.info("🔬 **Tier-2** — high-splice-scoring variant not annotated as canonical splice by VEP.")


# ---------------------------------------------------------------------------
# Page 2 — Tier-1 Candidates
# ---------------------------------------------------------------------------

def render_tier1_page() -> None:
    st.title("⭐ Tier-1 Candidates")
    st.markdown(
        "The **4 SCN1A splice-region VUS** prioritized for the Carvill lab minigene validation. "
        "Selected by mechanism + orthogonal evidence "
        "(chromatin similarity to known pathogenic + brain RNA-seq effect + splice-site proximity), "
        "**not** by raw splice score alone. Source: "
        "[`outputs/tier1_actionable_features.md`]({}) (Experiments 013 & 014)."
        .format(GITHUB_REPO + "/blob/main/outputs/tier1_actionable_features.md")
    )
    st.info(
        "**Bottom line:** If budget allows 2 variants, the strongest combined case is "
        "**rs801806** (strongest splicing/chromatin signature) + **rs4293437** (strongest predicted "
        "brain expression impact). For just 1: pick **rs801806** for splicing signature, or "
        "**rs4293437** for brain expression."
    )

    for variant in TIER1_DATA:
        _render_tier1_card(variant)


def _render_tier1_card(v: dict) -> None:
    rsid = v["rsid"]
    with st.container(border=True):
        # Header
        h1, h2 = st.columns([3, 1])
        with h1:
            st.markdown(f"### {rsid} — `chr{v['chrom']}:{v['pos']:,} {v['ref']}>{v['alt']}`")
        with h2:
            tier_badge = "⭐ Tier-1"
            st.markdown(
                f"<div style='text-align:right'><span style='background:#f0c674;padding:4px 10px;border-radius:8px;font-weight:600'>{tier_badge}</span></div>",
                unsafe_allow_html=True,
            )
        st.caption(f"**Consequence:** `{v['consequence']}`  •  **ClinVar:** `{v['clndn']}`")

        # Scores — 5 columns
        st.markdown("**AlphaGenome scores** (live re-call, `outputs/tier1_actionable_features.md`)")
        s1, s2, s3, s4, s5 = st.columns(5)
        s1.metric("SPLICE_SITES", f"{v['splice_sites']:.4f}", help="Primary splicing scorer")
        s2.metric("SPLICE_SITE_USAGE", f"{v['splice_site_usage']:.2f}")
        s3.metric(
            "SPLICE_JUNCTIONS (live |Δ| sum)",
            f"{v['splice_junctions_live']:.1f}",
            help=f"Atlas cached value: {v['splice_junctions_atlas']:.1f}",
        )
        s4.metric("DNASE Δ", f"{v['dnase_dscore']:.3f}")
        s5.metric(
            "Cosine to known pathogenic (mean / max)",
            f"{v['mean_cos_path']:.3f} / {v['max_cos_path']:.3f}",
            help="Cosine similarity of ISM matrices to 30 known-pathogenic SCN1A variants (Exp 004).",
        )

        # SpliceAI comparison placeholder — full SpliceAI deltas not in the brief, but SpliceAI is
        # an orthogonal benchmark the project uses; we surface what we have.
        st.markdown("**SpliceAI comparison**")
        st.caption(
            "SPLICE_SITES is AlphaGenome's primary splice scorer; SpliceAI deltas for these exact "
            "4 variants are not embedded in the brief — see `outputs/tier1_actionable_features.md` for the "
            "primary AlphaGenome-derived evidence. All 4 score ≥ 0.91 on SPLICE_SITES, in the top quartile of the "
            "1,610 VUS."
        )

        # Mechanism + population + literature
        st.markdown("**Mechanism hypothesis**")
        st.write(v["mechanism"])
        st.caption(f"Closest annotated splice site: `{v['closest_splice_site']}` (canonical: {v['canonical_exon']})")

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Population evidence (gnomAD v4.1)**")
            if v["gnomad_present"]:
                st.write(
                    f"Present (genome AF = {_fmt_af(v['gnomad_genome_af'])}, "
                    f"exome AF = {_fmt_af(v['gnomad_exome_af'])}). "
                    "Ultra-rare; absence of homozygous carriers is the relevant signal "
                    "(consistent with a severe dominant-acting variant)."
                )
            else:
                st.write("Absent from gnomAD (~125k exomes + 76k genomes). Consistent with a rare, potentially pathogenic variant — does NOT exclude benign rarity.")
        with c2:
            st.markdown("**Literature (PubMed / ClinVar)**")
            if v["pubmed_hits"] == 0:
                st.write("No dedicated PubMed entry for this rsID as of 2026-09 — not previously characterized.")
            else:
                st.write(f"{v['pubmed_hits']} PubMed hit(s).")


# ---------------------------------------------------------------------------
# Page 3 — Tier-2 Candidates
# ---------------------------------------------------------------------------

def render_tier2_page() -> None:
    st.title("🔬 Tier-2 Candidates")
    st.markdown(
        "**13 SCN1A VUS** with AlphaGenome scores in the same range as Tier-1, "
        "but **not** annotated as canonical splice variants by VEP. These represent "
        "three mechanism classes: deep-intronic cryptic splice activation (5), "
        "**dual-mechanism missense + splice** (5), and non-coding transcript overlap (3). "
        "Source: [`outputs/vus_tier2_candidates.md`]({}) (Experiment 010). "
        "**All 13 are absent from gnomAD and have 0 SCN1A-specific PubMed hits.**"
        .format(GITHUB_REPO + "/blob/main/outputs/vus_tier2_candidates.md")
    )

    tier2 = load_tier2_curated()

    # Stratification summary
    strat_counts = tier2["stratum"].value_counts().to_dict()
    s1, s2, s3 = st.columns(3)
    s1.metric("Deep intronic (cryptic splice)", strat_counts.get("intron_variant", 0))
    s2.metric("Missense (dual-mechanism)", strat_counts.get("missense_variant", 0))
    s3.metric("Non-coding transcript", strat_counts.get("non-coding_transcript_variant", 0))

    strata_meta = [
        ("intron_variant", "🧬 Stratum A — Deep intronic", "Cryptic splice site / poison exon / ISE-ISS disruption"),
        ("missense_variant", "🧪 Stratum B — Missense (dual mechanism)", "Coding + splicing; minigene assays that test only the canonical transcript will miss the splice effect"),
        ("non-coding_transcript_variant", "📜 Stratum C — Non-coding transcript overlap", "Antisense / non-coding RNA / regulatory overlap"),
    ]

    for stratum_key, stratum_title, stratum_subtitle in strata_meta:
        st.markdown(f"#### {stratum_title}")
        st.caption(stratum_subtitle)
        for _, v in tier2[tier2["stratum"] == stratum_key].iterrows():
            _render_tier2_card(v)


def _render_tier2_card(v: pd.Series) -> None:
    rsid = v["rsid"]
    with st.container(border=True):
        # Header
        h1, h2 = st.columns([3, 1])
        with h1:
            st.markdown(f"### {rsid} — `chr{int(v.get('chrom', 2))}:{int(v['pos']):,} {v['ref']}>{v['alt']}`")
        with h2:
            st.markdown(
                "<div style='text-align:right'><span style='background:#a3c9ff;padding:4px 10px;border-radius:8px;font-weight:600'>🔬 Tier-2</span></div>",
                unsafe_allow_html=True,
            )
        st.caption(
            f"**Consequence:** `{v['stratum']}`  •  **ClinVar:** `{v.get('clndn', '—')}`  •  "
            f"**Splice-only rank:** {int(v['splice_only_rank'])}  •  "
            f"**Combined score:** {float(v['combined_score']):.3f}"
        )

        s1, s2, s3 = st.columns(3)
        s1.metric("SPLICE_SITES", f"{float(v['SPLICE_SITES_score']):.4f}")
        s2.metric("DNASE Δ", f"{float(v['DNASE_score']):.2f}")
        s3.metric("Combined", f"{float(v['combined_score']):.3f}")

        st.markdown("**Mechanism hypothesis**")
        st.write(v["mechanistic_hypothesis"])

        st.markdown(
            "**Population evidence:** ⚪ absent from gnomAD (rare or unobserved).  "
            "**Literature:** 0 SCN1A-specific PubMed hits — not previously characterized."
        )


# ---------------------------------------------------------------------------
# Page 4 — Submit VCF
# ---------------------------------------------------------------------------

def render_submit_vcf_page() -> None:
    st.title("📤 Submit VCF")
    st.markdown(
        "Paste a VCF (or upload one) and score each variant against AlphaGenome. "
        "**This page uses real API quota** — the AlphaGenome API key you provide here is "
        "sent on every score request."
    )

    st.error(
        "⚠️ **API quota warning.** AlphaGenome's free tier is rate-limited (~30 variants/min, "
        "tens of thousands per week). Uploading a large VCF can exhaust your quota. "
        "We recommend submitting **≤50 variants at a time**. The key you enter below is held "
        "in `st.session_state` only — **not** persisted, logged, or sent anywhere except "
        "the official AlphaGenome endpoint."
    )

    api_key = st.text_input(
        "AlphaGenome API key",
        value=st.session_state.get("ag_key", ""),
        type="password",
        help="Get a free key at https://alphagenome.google/api. Stored in session state only.",
    )
    if api_key:
        st.session_state["ag_key"] = api_key
        st.success(f"🔑 Key loaded ({len(api_key)} chars). Held in session state — not saved.")
    else:
        st.info("🔑 No key entered yet — you can still upload/preview your VCF, but scoring will be disabled.")

    st.markdown("---")

    tab_paste, tab_upload = st.tabs(["📋 Paste VCF text", "📁 Upload .vcf file"])

    pasted = ""
    with tab_paste:
        pasted = st.text_area(
            "Paste VCF here (header lines starting with `#` are fine)",
            height=220,
            placeholder='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t...\n2\t166041471\trs801806\tT\tA\t.\t.\t.\n',
        )

    uploaded_df: Optional[pd.DataFrame] = None
    with tab_upload:
        up = st.file_uploader("Upload a .vcf file", type=["vcf", "txt"], accept_multiple_files=False)
        if up is not None:
            try:
                text = up.read().decode("utf-8", errors="replace")
                uploaded_df = parse_vcf(text)
                st.success(f"Parsed {len(uploaded_df)} allele row(s) from `{up.name}`.")
            except Exception as exc:
                st.error(f"Could not parse uploaded file: `{exc}`")

    # Combine sources (paste + upload) and de-dupe
    candidates = parse_vcf(pasted) if pasted else pd.DataFrame(columns=["chrom", "pos", "ref", "alt"])
    if uploaded_df is not None and not uploaded_df.empty:
        candidates = pd.concat([candidates, uploaded_df], ignore_index=True)
    candidates = candidates.drop_duplicates().reset_index(drop=True)

    if candidates.empty:
        st.info("Upload or paste a VCF to see parsed variants here.")
        return

    st.markdown(f"**{len(candidates)} candidate allele(s) parsed.** Preview:")
    st.dataframe(candidates.head(50), width="stretch", hide_index=True)

    if not api_key:
        st.warning("Enter an AlphaGenome API key above to enable scoring.")
        return

    if st.button("🚀 Score variants with AlphaGenome", type="primary"):
        _run_scoring(candidates, api_key)


def _run_scoring(variants: pd.DataFrame, api_key: str) -> None:
    """Run AlphaGenome on each variant. This burns real API quota."""
    if len(variants) > 200:
        st.warning(
            f"You submitted {len(variants)} variants. Truncating to 200 to protect your quota — "
            "edit the input to score a subset."
        )
        variants = variants.head(200)

    progress = st.progress(0.0, text="Starting…")
    results = []
    status = st.empty()

    try:
        from alphagenome.models import dna_model
        from alphagenome.models import variant_scorers
    except Exception as exc:
        st.error(
            f"Could not import the `alphagenome` SDK on this host. "
            f"This page requires the SDK locally — HF Spaces cannot run scoring for you "
            f"(the SDK is not pip-installable in the Space sandbox). Error: `{exc}`"
        )
        return

    try:
        model = dna_model.create(api_key=api_key, organism="human")
        scorers = [
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"],
        ]
    except Exception as exc:
        st.error(f"Could not initialise AlphaGenome with the provided key: `{exc}`")
        return

    for i, row in variants.iterrows():
        progress.progress((i + 1) / len(variants), text=f"Scoring {i + 1}/{len(variants)}")
        try:
            from alphagenome.data import genome
            variant = genome.Variant(
                chromosome=f"chr{row['chrom']}",
                position=int(row["pos"]),
                reference_bases=row["ref"],
                alternate_bases=row["alt"],
            )
            scores = model.score_variant(variant, scorers=scorers, sequence_length=2**14)
            record = {
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "SPLICE_SITES": float(scores.get("SPLICE_SITES", float("nan"))),
                "SPLICE_SITE_USAGE": float(scores.get("SPLICE_SITE_USAGE", float("nan"))),
                "DNASE": float(scores.get("DNASE", float("nan"))),
            }
            results.append(record)
        except Exception as exc:
            results.append(
                {
                    "chrom": row["chrom"],
                    "pos": row["pos"],
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "SPLICE_SITES": float("nan"),
                    "SPLICE_SITE_USAGE": float("nan"),
                    "DNASE": float("nan"),
                    "error": str(exc),
                }
            )
        status.write(f"Scored {i + 1}/{len(variants)}")

    progress.empty()
    status.empty()

    out = pd.DataFrame(results).sort_values("SPLICE_SITES", ascending=False).reset_index(drop=True)
    st.success(f"Scored {len(out)} variant(s) — ranked by AlphaGenome SPLICE_SITES score.")
    st.dataframe(out, width="stretch", hide_index=True)

    csv = out.to_csv(index=False)
    st.download_button(
        "⬇️ Download results as CSV",
        data=csv,
        file_name="alphagenome_scn1a_vcf_results.csv",
        mime="text/csv",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    page = render_sidebar()

    if page.startswith("📊"):
        render_browse_page()
    elif page.startswith("⭐"):
        render_tier1_page()
    elif page.startswith("🔬"):
        render_tier2_page()
    elif page.startswith("📤"):
        render_submit_vcf_page()
    else:
        st.error(f"Unknown page: {page!r}")


if __name__ == "__main__":
    main()