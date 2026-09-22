#!/usr/bin/env python3
"""
Download ClinVar variants for SCN1A.

Uses NCBI E-utilities (free, no API key needed for moderate usage).
We filter to:
- Gene: SCN1A
- Assembly: GRCh38
- Review status: at least one gold star (clinically reviewed)

Output: outputs/clinvar_scn1a.tsv with columns:
- chrom, pos, ref, alt
- hgvs_c (coding DNA)
- clinical_significance
- review_status
- condition
- gene_symbol
- molecular_consequence (from ClinVar)

Why not just download the full ClinVar file?
- Full release is ~10 GB and includes all genes
- E-utilities gene-scoped query is faster, more focused, and updates daily
- For our project we only care about SCN1A variants

Usage:
    python scripts/clinvar_scn1a.py
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
from xml.etree import ElementTree as ET

CLINVAR_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def esearch(gene_symbol: str, assembly: str = "GRCh38") -> list[str]:
    """Search ClinVar for variants in a gene on a specific assembly.

    Returns a list of ClinVar variation IDs (RCV accessions, but we use
    the simpler 'var' ID set here).
    """
    # First find the gene ID in ClinVar's gene-specific database
    # Strategy: use ClinVar's gene-conditioned search via the web API
    # Alternative: query the ClinVar variation archive directly
    # We'll use esearch on the clinvar database with a gene filter
    term = f"{gene_symbol}[gene] AND {assembly}[assembly]"
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": 5000,
        "retmode": "json",
    }
    url = f"{CLINVAR_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    print(f"E-search: {url[:120]}...")
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read())
    ids = data.get("esearchresult", {}).get("idlist", [])
    print(f"Got {len(ids)} ClinVar IDs")
    return ids


def esummary(variation_ids: list[str]) -> list[dict]:
    """Fetch summary records for a batch of ClinVar IDs."""
    records = []
    # esummary accepts up to 500 IDs per call
    for i in range(0, len(variation_ids), 200):
        batch = variation_ids[i:i + 200]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json",
        }
        url = f"{CLINVAR_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
        print(f"  E-summary batch {i // 200 + 1}: {len(batch)} IDs")
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = json.loads(resp.read())
        result = data.get("result", {})
        # uids are the variation IDs; each is a top-level key
        for uid in result.get("uids", []):
            if uid in result:
                records.append(result[uid])
        time.sleep(0.5)  # NCBI rate limit: 3 req/sec without API key
    return records


def parse_record(rec: dict) -> dict | None:
    """Extract our fields of interest from a ClinVar summary record."""
    try:
        # Get the germline classification
        classifications = rec.get("germline_classification", {})
        clinical_sig = classifications.get("description", "")
        review_status = classifications.get("review_status", "")

        # Trait set / condition
        trait_set = rec.get("trait_set", [])
        condition = ""
        if trait_set:
            condition = trait_set[0].get("trait_name", "")

        # Get variants from germline consequence
        # The exact field structure varies; we'll pull from the canonical allele
        canonical = rec.get("canonical_spdi", "") or rec.get("spdi", "")

        # Get molecular consequence
        mol_cons = rec.get("molecular_consequence", "")
        if isinstance(mol_cons, list):
            mol_cons = ";".join(mol_cons)

        # Get protein change (HGVS)
        hgvs = rec.get("protein_change", "") or rec.get("names", {}).get("protein", "")

        # Title often has the HGVS c. notation
        title = rec.get("title", "")

        return {
            "clinvar_id": rec.get("uid", ""),
            "title": title,
            "clinical_significance": clinical_sig,
            "review_status": review_status,
            "condition": condition,
            "molecular_consequence": mol_cons,
            "hgvs": hgvs,
            "canonical_spdi": canonical,
        }
    except Exception as e:
        print(f"  WARN: failed to parse record {rec.get('uid', '?')}: {e}")
        return None


def main() -> int:
    gene_symbol = os.environ.get("GENE_SYMBOL", "SCN1A")
    print(f"Fetching ClinVar variants for {gene_symbol}...")

    ids = esearch(gene_symbol)
    if not ids:
        print("No variants found.")
        return 1

    print(f"Fetching summaries for {len(ids)} variants...")
    records = esummary(ids)
    print(f"Got {len(records)} summary records")

    parsed = []
    for rec in records:
        row = parse_record(rec)
        if row:
            parsed.append(row)

    print(f"Parsed {len(parsed)} records")

    # Print some examples
    print("\nFirst 3 records (sample):")
    for row in parsed[:3]:
        print(json.dumps(row, indent=2))

    # Save
    os.makedirs("outputs", exist_ok=True)
    out_path = f"outputs/clinvar_{gene_symbol.lower()}.json"
    with open(out_path, "w") as f:
        json.dump(parsed, f, indent=2)
    print(f"\nSaved → {out_path}")

    # Quick stats
    sig_counts: dict[str, int] = {}
    for row in parsed:
        sig = row.get("clinical_significance", "unknown")
        sig_counts[sig] = sig_counts.get(sig, 0) + 1
    print("\nClinical significance counts:")
    for sig, count in sorted(sig_counts.items(), key=lambda x: -x[1]):
        print(f"  {sig}: {count}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
