#!/usr/bin/env python3
"""
Known SCN1A noncoding pathogenic variants — positive controls.

Curated from:
- Carvill et al. 2018 (AJHG) — original SCN1A poison exon paper, 7 variants
- Sparber et al. 2023 (Hum Genet) — minigene assay paper, 18 deep-intronic SNVs
- Frontiers 2025 review (Zhang et al.) — mentions c.4853-25 T>A as canonical example

Note: Many of these were originally reported in HGVS / transcript-relative
coordinates. We provide BOTH transcript coordinates (c. notation) and
genomic (chr:pos) coordinates where possible. Genomic coordinates need to
be verified against hg38 before use; some older papers use hg19.

If a variant only has transcript-relative coordinates, we mark it as
'needs_lift' and skip from genomic analysis (Atlas / API need chr:pos).
"""

from __future__ import annotations

# c. notation → expected mechanism, sourced from primary literature
KNOWN_PATHOGENIC_VARIANTS = [
    # === From Carvill et al. 2018 (AJHG) — 7 noncoding variants ===
    {
        "source": "Carvill 2018",
        "hgvs_c": "c.4853-25 T>A",
        "mechanism": "Activates cryptic splice site in intron 20",
        "location": "intron 20 (near poison exon 20N)",
        "patient_count": "multiple",
    },
    # The other 6 from Carvill 2018 are reported in the paper but specific
    # c. notation is scattered. Add them as we lift coordinates from the paper.
    {
        "source": "Carvill 2018",
        "hgvs_c": "(see paper Table S1)",
        "mechanism": "20N poison exon inclusion via SRSF1 motif disruption",
        "location": "intron 20 flanking 20N",
        "patient_count": "7 probands total",
    },
    # === From Sparber et al. 2023 (Hum Genet) ===
    # 18 deep-intronic SNVs tested in minigene assay; specific variants
    # listed in Table 1 of the paper. Add as we extract coordinates.
    {
        "source": "Sparber 2023",
        "hgvs_c": "(18 deep-intronic SNVs, see paper Table 1)",
        "mechanism": "Cryptic splice site activation / poison exon inclusion",
        "location": "introns 4, 6, 20, 21, 22, 24 (multiple poison exons)",
        "patient_count": "tested in vitro",
    },
    # === From Zhang 2025 Frontiers review — canonical example ===
    {
        "source": "Zhang 2025 (review)",
        "hgvs_c": "c.4853-25 T>A",
        "mechanism": "Partial exon skipping / intron retention → mild phenotype",
        "location": "intron 20 (near poison exon 20N)",
        "patient_count": "multiple",
    },
]

# Known SCN1A poison exons (from Carvill 2018 + Sparber 2023)
# These are the regions where regulatory variants are most likely pathogenic.
KNOWN_POISON_EXONS = [
    {"name": "20N", "intron": 20, "transcript": "NM_001165963.4"},
    # Additional poison exons identified by Sparber 2023
    # (specific names need to be lifted from the paper)
]


def main():
    print("SCN1A known pathogenic noncoding variants (literature)")
    print("=" * 60)
    for v in KNOWN_PATHOGENIC_VARIANTS:
        print(f"\n[{v['source']}]")
        print(f"  HGVS: {v['hgvs_c']}")
        print(f"  Mechanism: {v['mechanism']}")
        print(f"  Location: {v['location']}")
        print(f"  Patients: {v['patient_count']}")

    print("\n" + "=" * 60)
    print("Known SCN1A poison exons:")
    for pe in KNOWN_POISON_EXONS:
        print(f"  {pe['name']} (intron {pe['intron']}, transcript {pe['transcript']})")

    print("\n" + "=" * 60)
    print("NEXT STEP: extract specific c. notation + lift to genomic coordinates")
    print("from Carvill 2018 Table S1 and Sparber 2023 Table 1.")
    print("This requires PDF extraction from the papers.")


if __name__ == "__main__":
    main()
