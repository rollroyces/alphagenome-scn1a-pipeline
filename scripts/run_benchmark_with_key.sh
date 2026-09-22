#!/usr/bin/env bash
# Helper: run benchmark with key from secure prompt (not chat).
# Usage: bash scripts/run_benchmark_with_key.sh

set -e

KEY_FILE=".alphagenome_key"
chmod 600 "$KEY_FILE" 2>/dev/null || true

if [ ! -s "$KEY_FILE" ]; then
    echo "Enter your AlphaGenome API key (input hidden, not stored in chat):"
    read -s KEY
    echo "$KEY" > "$KEY_FILE"
    chmod 600 "$KEY_FILE"
    echo "Key saved to $KEY_FILE (mode 600, owned by you)."
    echo "To clear later: rm $KEY_FILE"
fi

export ALPHAGENOME_API_KEY="$(cat $KEY_FILE)"
echo "Key loaded (length: $(echo -n $ALPHAGENOME_API_KEY | wc -c | tr -d ' '))"

cd "$(dirname "$0")/.."
source .venv/bin/activate
python scripts/benchmark_scn1a_live_api.py
