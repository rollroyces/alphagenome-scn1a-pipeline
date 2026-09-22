#!/usr/bin/env python3
"""
Smoke test for the AlphaGenome API key.

Run this BEFORE opening the tutorial notebooks. It verifies:
1. The SDK imports cleanly
2. The API key is set
3. The model endpoint responds
4. A trivial prediction round-trips

Usage:
    export ALPHAGENOME_API_KEY=...
    python scripts/smoke_test.py

Expected output: a 1-row DataFrame with the model's prediction, then "OK".
Expected failure modes:
    - Missing key → clear error
    - Wrong key → clear error from the API
    - No network → timeout error
"""

from __future__ import annotations

import os
import sys

import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY environment variable is not set.")
        print("Get a free key at https://alphagenome.google/api (non-commercial use).")
        return 1

    print(f"API key length: {len(api_key)} chars (looks plausible)")
    print("Creating dna_model client...")

    try:
        dna_model = dna_client.create(api_key)
    except Exception as e:
        print(f"ERROR: Failed to create client: {type(e).__name__}: {e}")
        return 2

    print(f"Client created. Sequence lengths supported: {list(dna_client.SUPPORTED_SEQUENCE_LENGTHS.keys())}")
    print(f"Output types: {[o.name for o in dna_client.OutputType]}")

    print("\nMaking trivial prediction: 1bp of 'A' padded to 16Kb, DNase in lung...")
    sequence = "A".center(dna_client.SEQUENCE_LENGTH_16KB, "N")

    try:
        output = dna_model.predict_sequence(
            sequence=sequence,
            requested_outputs=[dna_client.OutputType.DNASE],
            ontology_terms=["UBERON:0002048"],  # Lung
        )
    except Exception as e:
        print(f"ERROR: predict_sequence failed: {type(e).__name__}: {e}")
        return 3

    dnase = output.dnase
    print(f"\nGot TrackData:")
    print(f"  values shape: {dnase.values.shape}")
    print(f"  dtype: {dnase.values.dtype}")
    print(f"  metadata columns: {list(dnase.metadata.columns)}")
    print(f"\nFirst few metadata rows:")
    print(dnase.metadata.head().to_string())
    print(f"\nFirst 5 values (raw prediction, any number is fine — 'A' isn't a real sequence):")
    print(dnase.values[:5])

    print("\nOK — AlphaGenome API is reachable and your key works.")
    print("Next: open tutorials/alphagenome-upstream/colabs/quick_start.ipynb in Jupyter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
