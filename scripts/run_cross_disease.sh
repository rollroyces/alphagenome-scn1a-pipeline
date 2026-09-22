#!/usr/bin/env bash
# Run cross-disease benchmark for all 5 genes sequentially.
# Total runtime: ~20-25 minutes.
set -e
cd /Users/hermes/projects/alphagenome-work
source .venv/bin/activate
export ALPHAGENOME_API_KEY=$(cat .alphagenome_key)

for gene in scn1a scn2a mecp2 cftr dmd; do
    echo ""
    echo "=========================================="
    echo "Scoring $gene..."
    echo "=========================================="
    python scripts/score_gene_benchmark.py \
        --input outputs/clinvar_${gene}_benchmark.tsv \
        --output outputs/cross_disease_${gene}_raw.csv
done

echo ""
echo "=========================================="
echo "All genes scored. Aggregating..."
echo "=========================================="
python scripts/aggregate_cross_disease.py
