#!/usr/bin/env python3
"""
Experiment 008 — KCNQ2 extraction + stratification.

Extract SNV variants for KCNQ2 from local ClinVar VCF, stratify by clnsig,
dump to a TSV that we'll score with the live API.

Output: outputs/_kcnq2_stratified.tsv (one row per variant, columns:
chrom, pos, ref, alt, rsid, clnsig, clnsig_category, gene_label)
"""

from __future__ import annotations

import csv
import sys

import pysam

# KCNQ2 — MANE Select ENST00000356457 / NM_172107.4 (hg38, plus strand)
GENE = "KCNQ2"
GENE_ID = "3785"
CHROM = "20"
START = 63400679
END = 63472909
VCF = "data/clinvar_grch38.vcf.gz"
TABIX = "data/clinvar_grch38.vcf.gz.tbi"
OUT_TSV = "outputs/_kcnq2_stratified.tsv"


def parse_info(s: str) -> dict[str, str]:
    out = {}
    for entry in s.split(";"):
        if "=" in entry:
            k, v = entry.split("=", 1)
            out[k] = v
    return out


def get_gene_names(info: dict) -> list[str]:
    g = info.get("GENEINFO", "")
    out = []
    for entry in g.split("|"):
        if ":" in entry:
            out.append(entry.split(":", 1)[0])
        elif entry:
            out.append(entry)
    return out


def clnsig_category(clnsig: str) -> str:
    s = clnsig.lower()
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        # includes "likely_pathogenic"
        return "pathogenic"
    if "benign" in s and "pathogenic" not in s:
        # includes "likely_benign"
        return "benign"
    if "uncertain" in s or "vus" in s:
        return "uncertain"
    return "other"


def main() -> int:
    tb = pysam.TabixFile(VCF, index=TABIX)
    region = f"{CHROM}:{START}-{END}"
    recs = list(tb.fetch(region))
    tb.close()
    print(f"Total records in {region}: {len(recs)}")

    rows = []
    for line in recs:
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        info = parse_info(fields[7])
        if GENE not in get_gene_names(info):
            continue
        ref = fields[3]
        alts = fields[4].split(",") if fields[4] != "." else []
        # SNV only
        if len(ref) != 1:
            continue
        # Take first SNV alt (skip indels)
        snv_alt = None
        for a in alts:
            if len(a) == 1:
                snv_alt = a
                break
        if snv_alt is None:
            continue

        rsid = fields[2] if fields[2] != "." else ""
        clnsig = info.get("CLNSIG", "")
        cat = clnsig_category(clnsig)
        rows.append({
            "chrom": CHROM,
            "pos": int(fields[1]),
            "ref": ref,
            "alt": snv_alt,
            "rsid": rsid,
            "clnsig": clnsig,
            "clnsig_category": cat,
        })

    print(f"KCNQ2 SNVs: {len(rows)}")

    # Counts by category
    from collections import Counter
    cnt = Counter(r["clnsig_category"] for r in rows)
    print("By clnsig_category:")
    for k, v in sorted(cnt.items(), key=lambda kv: -kv[1]):
        print(f"  {k}: {v}")

    # Write TSV
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "rsid", "clnsig", "clnsig_category"])
        for r in rows:
            w.writerow([r["chrom"], r["pos"], r["ref"], r["alt"],
                        r["rsid"], r["clnsig"], r["clnsig_category"]])
    print(f"Saved → {OUT_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
