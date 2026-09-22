#!/usr/bin/env python3
"""
Converted-from-notebook version of quick_start.ipynb for fast local iteration.

Why this exists: Jupyter is great for exploration but slow when you just want
to verify a pipeline works. This script runs the same code as the first half
of quick_start.ipynb and writes a self-contained report.

It does NOT replace the full notebook — it gives you a fast feedback loop
for "does the API work, do I understand the output shapes, can I do batch calls."

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/quickstart_script.py
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")  # headless backend, no display required
import matplotlib.pyplot as plt
import pandas as pd

from alphagenome.data import gene_annotation, genome, transcript as transcript_utils
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components


def section(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Set ALPHAGENOME_API_KEY first.")
        return 1

    dna_model = dna_client.create(api_key)

    section("Step 1: Trivial prediction — padded sequence, DNase in lung")
    output = dna_model.predict_sequence(
        sequence="GATTACA".center(dna_client.SEQUENCE_LENGTH_1MB, "N"),
        requested_outputs=[dna_client.OutputType.DNASE],
        ontology_terms=["UBERON:0002048"],  # Lung
    )
    dnase = output.dnase
    print(f"DNase values shape: {dnase.values.shape}")
    print(f"  → (sequence_length, num_tracks) = (1,048,576, 1)")
    print(f"  → values are float32, range typically [0, ~5]")
    print(f"  → first 3 values: {dnase.values[:3].flatten()}")
    print(f"\nMetadata columns: {list(dnase.metadata.columns)}")

    section("Step 2: Multi-tissue prediction")
    output = dna_model.predict_sequence(
        sequence="GATTACA".center(dna_client.SEQUENCE_LENGTH_1MB, "N"),
        requested_outputs=[dna_client.OutputType.CAGE, dna_client.OutputType.DNASE],
        ontology_terms=[
            "UBERON:0002048",  # Lung
            "UBERON:0000955",  # Brain
        ],
    )
    print(f"DNASE shape: {output.dnase.values.shape}  → (1Mb, 2 tissues)")
    print(f"CAGE shape:  {output.cage.values.shape}  → (1Mb, 4 = 2 tissues × 2 strands)")
    print(f"\nCAGE metadata (each row = one track):")
    print(output.cage.metadata.to_string())

    section("Step 3: Reference genome interval — RNA-seq at CYP2B6")
    print("Loading GENCODE v46 GTF (~150 MB, takes ~30s first time)...")
    gtf = pd.read_feather(
        "https://storage.googleapis.com/alphagenome/reference/gencode/"
        "hg38/gencode.v46.annotation.gtf.gz.feather"
    )
    print(f"GTF loaded: {len(gtf):,} annotation rows")

    gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
    gtf_transcripts = gene_annotation.filter_to_mane_select_transcript(gtf_transcripts)
    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)
    print(f"MANE Select transcripts: {len(gtf_transcripts):,}")

    interval = gene_annotation.get_gene_interval(gtf, gene_symbol="CYP2B6")
    print(f"CYP2B6 interval: {interval}")
    interval = interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    print(f"Resized to {interval.width:,} bp (1 Mb centered on gene)")

    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
        ontology_terms=["UBERON:0001114"],  # Right liver lobe
    )
    print(f"RNA-seq shape: {output.rna_seq.values.shape}")
    print(f"Metadata:")
    print(output.rna_seq.metadata.to_string())

    transcripts = transcript_extractor.extract(interval)
    print(f"Transcripts in interval: {len(transcripts)}")

    fig, ax = plt.subplots(figsize=(12, 4))
    plot_components.plot(
        components=[
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(output.rna_seq),
        ],
        interval=output.rna_seq.interval,
        ax=ax,
    )
    plt.tight_layout()
    out_path = "figures/quickstart_cyp2b6_rnaseq.png"
    os.makedirs("figures", exist_ok=True)
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved locus plot → {out_path}")

    section("Step 4: Variant effect — chr22:36201698 A>C, RNA-seq in colon")
    variant = genome.Variant(
        chromosome="chr22",
        position=36201698,
        reference_bases="A",
        alternate_bases="C",
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    print(f"Variant: {variant}")
    print(f"Reference interval: {interval}")

    variant_output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
        ontology_terms=["UBERON:0001157"],  # Colon - Transverse
    )
    print(f"\nGot ref and alt predictions.")
    print(f"Ref RNA-seq shape: {variant_output.reference.rna_seq.values.shape}")
    print(f"Alt RNA-seq shape: {variant_output.alternate.rna_seq.values.shape}")

    transcripts = transcript_extractor.extract(interval)
    fig, ax = plt.subplots(figsize=(12, 5))
    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.OverlaidTracks(
                tdata={
                    "REF": variant_output.reference.rna_seq,
                    "ALT": variant_output.alternate.rna_seq,
                },
                colors={"REF": "dimgrey", "ALT": "red"},
            ),
        ],
        interval=variant_output.reference.rna_seq.interval.resize(2**15),
        annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
        ax=ax,
    )
    plt.tight_layout()
    out_path = "figures/quickstart_variant_ref_vs_alt.png"
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved ref/alt comparison → {out_path}")

    section("Step 5: Variant scoring — collapse ref/alt deltas to a per-gene score")
    variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
    print(f"Using scorer: {variant_scorer}")

    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=[variant_scorer],
    )
    print(f"Got {len(variant_scores)} AnnData objects (one per scorer)")

    adata = variant_scores[0]
    print(f"\nAnnData shape: {adata.X.shape}")
    print(f"  → (num_variants x num_tracks, but here num_variants=1)")
    print(f"\n.obs (variant metadata):")
    print(adata.obs.head().to_string())
    print(f"\n.var (track metadata, first 5):")
    print(adata.var.head().to_string())

    tidy = variant_scorers.tidy_scores([adata], match_gene_strand=True)
    print(f"\nTidy scores (long format, gene-matched):")
    print(tidy.head(10).to_string())

    section("DONE")
    print("If you reached here with no errors, you have:")
    print("  ✓ Working AlphaGenome API access")
    print("  ✓ Understood predict / score_variant loop")
    print("  ✓ Seen all key output types and their shapes")
    print("  ✓ Made your first locus plots")
    print("\nNext: open the visualization_modality_tour notebook to see all 11 modalities.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
