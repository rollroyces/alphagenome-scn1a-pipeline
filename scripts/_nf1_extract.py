#!/usr/bin/env python3
"""
Experiment 016 — NF1 extraction + stratification.

Extract SNV variants for NF1 from local ClinVar VCF, stratify by clnsig,
dump to a TSV that we'll score with the live API.

NF1 — GENCODE v46 (hg38):
  - chr17: 31,094,926 - 31,382,116 (gene span; ENSG00000196712.20)
  - Strand: plus
  - MANE Select transcript: ENST00000358273.9 (NF1-202, NM_000267.4)
  - NCBI Gene ID: 4763
  - Largest gene in the cross-disease benchmark by exon count (58 exons),
    2nd-largest by span (~287 kb), behind FBN1 (~237 kb but only 65 exons
    for FBN1's MANE Select).

Note on extraction:
  The pysam TabixFile.fetch() iterator auto-decodes bytes as ASCII, which
  fails on ClinVar records whose INFO carries non-ASCII text (e.g. the
  "Café-au-lait_macules" disease name in the NF1 region). We shell out
  to `tabix` and decode bytes ourselves in UTF-8. Output schema is
  identical to the other extract scripts (chrom, pos, ref, alt, rsid,
  clnsig, clnsig_category, mc).

Output: outputs/_nf1_stratified.tsv
"""

from __future__ import annotations

import csv
import subprocess
import sys

# NF1 — GENCODE v46 (hg38), ENSG00000196712
GENE = "NF1"
GENE_ID = "4763"
CHROM = "17"
START = 31_094_926
END = 31_382_116
VCF = "data/clinvar_grch38.vcf.gz"
TABIX = "data/clinvar_grch38.vcf.gz.tbi"
OUT_TSV = "outputs/_nf1_stratified.tsv"


def fetch_records(chrom: str, start: int, end: int) -> list[str]:
    """Fetch ClinVar records via the `tabix` CLI; decode bytes as UTF-8.

    Robust to non-ASCII entries that crash pysam's ASCII-only iterator.
    """
    out = subprocess.run(
        ["tabix", VCF, f"{chrom}:{start}-{end}"],
        capture_output=True,
        check=True,
    )
    text = out.stdout.decode("utf-8", errors="replace")
    return [ln for ln in text.split("\n") if ln]


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
    print(f"Fetching {GENE} region chr{CHROM}:{START}-{END} via tabix...")
    raw_lines = fetch_records(CHROM, START, END)
    print(f"Total records in region: {len(raw_lines)}")

    rows = []
    for line in raw_lines:
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
        mc = info.get("MC", "")
        rows.append({
            "chrom": CHROM,
            "pos": int(fields[1]),
            "ref": ref,
            "alt": snv_alt,
            "rsid": rsid,
            "clnsig": clnsig,
            "clnsig_category": cat,
            "mc": mc,
        })

    print(f"NF1 SNVs: {len(rows)}")

    # Counts by category
    from collections import Counter
    cnt = Counter(r["clnsig_category"] for r in rows)
    print("By clnsig_category:")
    for k, v in sorted(cnt.items(), key=lambda kv: -kv[1]):
        print(f"  {k}: {v}")

    # Counts of MC within pathogenic (so we can confirm the splice burden)
    mc_cnt = Counter()
    for r in rows:
        if r["clnsig_category"] != "pathogenic":
            continue
        for tag in [t.strip().lower() for t in r["mc"].split(",")]:
            if tag:
                mc_cnt[tag] += 1
    print("\nMolecular-consequence tags within pathogenic (top 15):")
    for k, v in mc_cnt.most_common(15):
        print(f"  {k}: {v}")

    # Write TSV
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "rsid",
                    "clnsig", "clnsig_category", "mc"])
        for r in rows:
            w.writerow([r["chrom"], r["pos"], r["ref"], r["alt"],
                        r["rsid"], r["clnsig"], r["clnsig_category"],
                        r["mc"]])
    print(f"\nSaved → {OUT_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
