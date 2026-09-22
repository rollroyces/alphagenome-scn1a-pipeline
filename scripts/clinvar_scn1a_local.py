#!/usr/bin/env python3
"""
Pull ClinVar variants for SCN1A from a local tabix-indexed ClinVar VCF.

This is the production version of the data acquisition step. Reads from
data/clinvar_grch38.vcf.gz (downloaded from NCBI), filters to SCN1A,
and produces a clean TSV ready for downstream analysis.

Output: outputs/clinvar_scn1a.tsv
- Columns: chrom, pos, ref, alt, rsid, clnsig, clndn, molecular_consequence, gene

ClinVar's CLNSIG uses descriptive strings, not codes. The mapping we use:
- Pathogenic / Likely_pathogenic / Pathogenic,_Low_Penetrance → pathogenic
- Benign / Likely_benign → benign
- Uncertain_significance / VUS → uncertain
- Conflicting_classifications → conflicting
- everything else → other

Usage:
    python scripts/clinvar_scn1a_local.py
"""

from __future__ import annotations

import csv
import os
import sys

import pysam


# SCN1A gene interval (hg38)
SCN1A_INTERVAL = ("2", 165_984_640, 166_182_806)
SCN1A_FLANK = 500_000  # regulatory variants can be far from gene body
SCN1A_GENE_ID = "6323"  # NCBI Gene ID for SCN1A


# ClinVar consequence ontology terms relevant to splicing / regulatory variants
SPLICING_RELATED_CONSEQUENCES = {
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_region_variant",
    "intron_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant",
}

PATHOGENIC_LABELS = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Pathogenic,_low_penetrance",
    "Pathogenic,_protective",
}

BENIGN_LABELS = {
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign",
}

UNCERTAIN_LABELS = {
    "Uncertain_significance",
    "Uncertain_significance,_conflicting_interpretations",
    "VUS-high",
}

CONFLICTING_LABELS = {
    "Conflicting_classifications_of_pathogenicity",
}


def parse_clinvar_record(rec_str: str) -> dict | None:
    """Parse a single tabix-returned record string into a structured dict."""
    if rec_str.startswith("#"):
        return None
    fields = rec_str.split("\t")
    if len(fields) < 8:
        return None
    chrom, pos, rsid, ref, alts, qual, filt, info = fields[:8]

    info_dict: dict[str, list[str]] = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info_dict[k] = v.split(",")
        else:
            info_dict[kv] = ["true"]

    return {
        "chrom": chrom,
        "pos": int(pos),
        "rsid": rsid,
        "ref": ref,
        "alts": alts.split(","),
        "info": info_dict,
    }


def normalize_clnsig(clnsig: str) -> str:
    """Collapse ClinVar's many descriptive labels into broad categories."""
    if clnsig in PATHOGENIC_LABELS:
        return "pathogenic"
    if clnsig in BENIGN_LABELS:
        return "benign"
    if clnsig in UNCERTAIN_LABELS:
        return "uncertain"
    if clnsig in CONFLICTING_LABELS:
        return "conflicting"
    return "other"


def parse_consequence(mc_str: str) -> str:
    """Extract the SO term from 'SO:0001623|intron_variant' format."""
    if "|" in mc_str:
        return mc_str.split("|", 1)[1]
    return mc_str


def is_scn1a(geneinfo: str) -> bool:
    """Check if SCN1A is in the GENEINFO field.

    ClinVar's GENEINFO uses pipe (not comma) to separate multiple genes:
        GENEINFO=SCN1A:6323|LOC102724058:102724058
    """
    if not geneinfo:
        return False
    for entry in geneinfo.split("|"):
        if ":" in entry:
            gene_name, gene_id = entry.split(":", 1)
            if gene_id == SCN1A_GENE_ID:
                return True
    return False


def main() -> int:
    vcf_path = "data/clinvar_grch38.vcf.gz"
    tbi_path = "data/clinvar_grch38.vcf.gz.tbi"
    if not os.path.exists(vcf_path) or not os.path.exists(tbi_path):
        print(f"ERROR: {vcf_path} or {tbi_path} not found.")
        print("Run: curl -o data/clinvar_grch38.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260913.vcf.gz")
        print("Run: curl -o data/clinvar_grch38.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260913.vcf.gz.tbi")
        return 1

    chrom, start, end = SCN1A_INTERVAL
    start_flank = max(0, start - SCN1A_FLANK)
    end_flank = end + SCN1A_FLANK
    region = f"{chrom}:{start_flank}-{end_flank}"

    print(f"Querying ClinVar for region {region} (SCN1A ± {SCN1A_FLANK:,} bp)...")
    tb = pysam.TabixFile(vcf_path, index=tbi_path)
    # tabix returns an iterator — materialize once, then reuse
    raw_records = list(tb.fetch(region))
    print(f"Got {len(raw_records):,} raw records in region")

    # Parse and filter to SCN1A
    scn1a_records = []
    all_records = []
    for rec_str in raw_records:
        parsed = parse_clinvar_record(rec_str)
        if parsed is None:
            continue
        all_records.append(parsed)
        geneinfo = parsed["info"].get("GENEINFO", [""])[0]
        if is_scn1a(geneinfo):
            scn1a_records.append(parsed)

    print(f"SCN1A-specific records (gene_id={SCN1A_GENE_ID}): {len(scn1a_records):,}")

    # Build output rows
    output_rows = []
    for rec in scn1a_records:
        info = rec["info"]
        clnsig = info.get("CLNSIG", [""])[0]
        clndn = info.get("CLNDN", [""])[0].replace("|", ",")
        mc_raw = info.get("MC", [""])[0]
        consequence = parse_consequence(mc_raw)
        # Take the first alt (multi-allelic sites handled separately if needed)
        alt = rec["alts"][0] if rec["alts"] else ""

        output_rows.append({
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": alt,
            "rsid": rec["rsid"],
            "clnsig_raw": clnsig,
            "clnsig_category": normalize_clnsig(clnsig),
            "clndn": clndn,
            "molecular_consequence": consequence,
            "is_splicing_related": "yes" if consequence in SPLICING_RELATED_CONSEQUENCES else "no",
            "review_status": info.get("CLNREVSTAT", [""])[0],
        })

    # Write TSV
    os.makedirs("outputs", exist_ok=True)
    out_path = "outputs/clinvar_scn1a.tsv"
    fieldnames = [
        "chrom", "pos", "ref", "alt", "rsid",
        "clnsig_raw", "clnsig_category",
        "clndn", "molecular_consequence",
        "is_splicing_related", "review_status",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)
    print(f"Saved → {out_path}")

    # Summary stats
    print("\n=== Summary ===")
    print(f"Total SCN1A records: {len(output_rows):,}")

    sig_counts: dict[str, int] = {}
    for row in output_rows:
        sig_counts[row["clnsig_category"]] = sig_counts.get(row["clnsig_category"], 0) + 1
    print("\nClinical significance:")
    for sig, count in sorted(sig_counts.items(), key=lambda x: -x[1]):
        print(f"  {sig}: {count}")

    splic_counts: dict[str, int] = {}
    for row in output_rows:
        splic_counts[row["molecular_consequence"]] = splic_counts.get(row["molecular_consequence"], 0) + 1
    print("\nMolecular consequence (top 15):")
    for mc, count in sorted(splic_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {mc}: {count}")

    # Specifically: pathogenic splicing-related variants (our positive controls)
    pathogenic_splicing = [
        r for r in output_rows
        if r["clnsig_category"] == "pathogenic" and r["is_splicing_related"] == "yes"
    ]
    print(f"\n★ Pathogenic + splicing-related variants: {len(pathogenic_splicing)} (gold-standard positive controls)")
    for r in pathogenic_splicing[:10]:
        print(f"  {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} | {r['molecular_consequence']} | {r['clndn'][:60]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
