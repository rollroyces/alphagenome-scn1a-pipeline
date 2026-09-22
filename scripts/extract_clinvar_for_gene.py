#!/usr/bin/env python3
"""
Gene-agnostic ClinVar benchmark extractor for AlphaGenome.

For a given gene, extract:
- Positive controls: pathogenic + splicing-related + SNV
- Negative controls: benign + intronic + SNV

Outputs: TSV with one row per variant.

Usage:
    python extract_clinvar_for_gene.py --gene SCN1A --chrom 2 --start 165984640 --end 166182806 --n-pos 200 --n-neg 350
"""

from __future__ import annotations

import argparse
import csv
import sys

import pysam


# Molecular consequences that count as "splicing-related"
# Includes both ClinVar MC labels and Sequence Ontology (SO) term names.
SPLICING_CONSEQUENCES = {
    # ClinVar labels
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
    # SO terms (numeric IDs and names)
    "SO:0001575",  # splice_donor_variant
    "SO:0001574",  # splice_acceptor_variant
    "splice donor variant",
    "splice acceptor variant",
    "splice region variant",
}

# Clinical significance categories
PATHOGENIC_LABELS = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
    "likely pathogenic",
}
BENIGN_LABELS = {
    "benign",
    "likely_benign",
    "benign/likely_benign",
    "likely benign",
}
UNCERTAIN_LABELS = {
    "uncertain",
    "uncertain_significance",
    "uncertain significance",
    "vus",
}


def parse_info(info_str: str) -> dict[str, list[str]]:
    """Parse a VCF INFO field into {key: [values]}."""
    out = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            k, v = entry.split("=", 1)
            out.setdefault(k, []).append(v)
        else:
            out.setdefault(entry, [""])
    return out


def parse_clinvar_record(line: str) -> dict:
    """Parse a tabix-returned ClinVar record string."""
    fields = line.split("\t")
    if len(fields) < 8:
        return {}
    info = parse_info(fields[7])
    return {
        "chrom": fields[0],
        "pos": int(fields[1]),
        "rsid": fields[2] if fields[2] != "." else "",
        "ref": fields[3],
        "alts": fields[4].split(",") if fields[4] != "." else [],
        "info": info,
    }


def get_consequence_list(info: dict) -> list[str]:
    """Get MC field as a list of consequence strings (handles pipe-delimited inner)."""
    mc = info.get("MC", [""])[0]
    return mc.split("|")


def get_gene_names(info: dict) -> list[str]:
    """Get the list of gene symbols in GENEINFO."""
    g = info.get("GENEINFO", [""])[0]
    out = []
    for entry in g.split("|"):
        if ":" in entry:
            sym = entry.split(":", 1)[0]
            out.append(sym)
        elif entry:
            out.append(entry)
    return out


def get_clnsig_category(info: dict) -> str:
    """Get a coarse category from CLNSIG."""
    raw = str(info.get("CLNSIG", [""])[0]).lower()
    if "conflicting" in raw:
        return "conflicting"
    if "pathogenic" in raw and "benign" not in raw:
        return "pathogenic"
    if "benign" in raw and "pathogenic" not in raw:
        return "benign"
    if "uncertain" in raw or "vus" in raw:
        return "uncertain"
    return "other"


def fetch_variants(chrom: str, start: int, end: int, gene_name: str,
                   vcf_path: str = "data/clinvar_grch38.vcf.gz",
                   tabix_path: str = "data/clinvar_grch38.vcf.gz.tbi") -> list[dict]:
    """Fetch all ClinVar variants for a gene, return list of dicts."""
    tb = pysam.TabixFile(vcf_path, index=tabix_path)
    # Try chr-prefixed first, fall back to no-prefix
    records = []
    for region in [f"chr{chrom}:{start}-{end}", f"{chrom}:{start}-{end}"]:
        try:
            records = list(tb.fetch(region))
            if records:
                break
        except ValueError:
            continue
    tb.close()

    out = []
    for line in records:
        rec = parse_clinvar_record(line)
        if not rec:
            continue
        if gene_name not in get_gene_names(rec["info"]):
            continue
        # Only SNVs (skip indels for AlphaGenome compatibility)
        if len(rec["ref"]) != 1 or any(len(a) != 1 for a in rec["alts"]):
            continue
        rec["clnsig_category"] = get_clnsig_category(rec["info"])
        rec["consequence_list"] = get_consequence_list(rec["info"])
        rec["molecular_consequence"] = ",".join(rec["consequence_list"])
        out.append(rec)
    return out


def make_benchmark(chrom: str, start: int, end: int, gene_name: str,
                   n_pos: int = 200, n_neg: int = 350, seed: int = 42) -> list[dict]:
    """Build a balanced benchmark from ClinVar for a gene."""
    import random
    random.seed(seed)

    variants = fetch_variants(chrom, start, end, gene_name)
    print(f"  Total SNV variants for {gene_name}: {len(variants)}")

    # Positives: pathogenic + has at least one splicing-related consequence
    positives = []
    seen = set()
    for v in variants:
        if v["clnsig_category"] != "pathogenic":
            continue
        has_spl = any(c in SPLICING_CONSEQUENCES for c in v["consequence_list"])
        if not has_spl:
            continue
        key = (v["chrom"], v["pos"], v["ref"], v["alts"][0])
        if key in seen:
            continue
        seen.add(key)
        positives.append({"variant": v, "label": "positive"})

    random.shuffle(positives)
    positives = positives[:n_pos]

    # Negatives: benign + intronic (with or without splice annotation)
    negatives = []
    seen = set()
    for v in variants:
        if v["clnsig_category"] != "benign":
            continue
        has_intron = any("intron" in c for c in v["consequence_list"])
        if not has_intron:
            continue
        key = (v["chrom"], v["pos"], v["ref"], v["alts"][0])
        if key in seen:
            continue
        seen.add(key)
        negatives.append({"variant": v, "label": "negative"})

    random.shuffle(negatives)
    negatives = negatives[:n_neg]

    benchmark = positives + negatives
    # Shuffle order to avoid API call order bias
    random.shuffle(benchmark)

    print(f"  Positives: {len(positives)}, Negatives: {len(negatives)}")
    return benchmark


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--chrom", required=True, help="Without 'chr' prefix")
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--n-pos", type=int, default=200)
    parser.add_argument("--n-neg", type=int, default=350)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    if args.out is None:
        args.out = f"outputs/clinvar_{args.gene.lower()}_benchmark.tsv"

    benchmark = make_benchmark(
        args.chrom, args.start, args.end, args.gene,
        n_pos=args.n_pos, n_neg=args.n_neg,
    )
    if len(benchmark) == 0:
        print("ERROR: No variants found")
        return 1
    labels = [b["label"] for b in benchmark]
    n_pos = sum(1 for l in labels if l == "positive")
    n_neg = sum(1 for l in labels if l == "negative")

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "chrom", "pos", "rsid", "ref", "alt",
                    "label", "clnsig_category", "molecular_consequence"])
        for b in benchmark:
            v = b["variant"]
            w.writerow([args.gene, v["chrom"], v["pos"], v["rsid"],
                        v["ref"], v["alts"][0], b["label"],
                        v["clnsig_category"], v["molecular_consequence"]])

    print(f"Saved → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
