#!/usr/bin/env bash
# Generic launcher: read key from .alphagenome_key, export, run script with args.
# Usage: bash scripts/_run_with_key.sh <script.py> [args...]
set -e
KEY=$(cat /Users/hermes/projects/alphagenome-work/.alphagenome_key)
export ALPHAGENOME_API_KEY="$KEY"
cd /Users/hermes/projects/alphagenome-work
source .venv/bin/activate
python "$@"
