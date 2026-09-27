# Makefile for alphagenome-scn1a-pipeline
# Reproducibility harness — re-verify existing outputs without burning the
# 4,000+ API calls needed to redo the full pipeline.

# ----- Config ----------------------------------------------------------------

PYTHON       ?= python3
VENV         ?= .venv
VENV_PY      := $(VENV)/bin/python
KEY_FILE     := .alphagenome_key
KEY_LAUNCHER := bash scripts/_run_with_key.sh

# Patterns of expected output files. Kept short; we only list outputs the
# project README / paper actually cites. Expand as new artefacts land.
EXPECTED_OUTPUTS := \
	outputs/benchmark_scn1a_live_api_raw.csv \
	outputs/benchmark_scn1a_results.csv \
	outputs/clinvar_scn1a.json \
	outputs/clinvar_scn1a.tsv \
	outputs/clinvar_scn1a_benchmark.tsv \
	outputs/vus_rescored.csv \
	outputs/vus_rescored_with_dnase.csv \
	outputs/vus_top_candidates.csv \
	outputs/vus_top_by_dnase.csv \
	outputs/vus_top_combined.csv \
	outputs/vus_high_impact_candidates.csv \
	outputs/vus_high_impact_with_dnase.csv \
	outputs/vus_high_impact_with_gnomad.csv \
	outputs/vus_tier2_candidates.csv \
	outputs/tier1_actionable_features.csv \
	outputs/tier1_actionable_features.md \
	outputs/cross_disease_metrics.csv \
	outputs/cross_disease_scn1a_raw.csv \
	outputs/cross_disease_summary.md

# Per-experiment artefacts we know exist.
EXP001_OUTPUTS := \
	research_notebook/experiments/001_ism_scn1a/metrics.csv
EXP011_OUTPUTS := \
	research_notebook/experiments/011_ism_tier1/metrics.csv
EXP013_OUTPUTS := \
	research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.csv \
	research_notebook/experiments/013_tier1_actionable/cryptic_splice_scores.json \
	outputs/tier1_actionable_features.csv

# ----- Targets ---------------------------------------------------------------

.PHONY: help check key smoke exp001 exp011 exp013 reproduce clean venv test

help:  ## Show this help.
	@awk 'BEGIN {FS = ":[^#]*## "} /^[a-zA-Z0-9_-]+:[^#]*## / {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

venv: ## Create the .venv if it doesn't exist (uv-managed, Python 3.13).
	@test -d $(VENV) || (uv venv --python 3.13 $(VENV) && uv pip install --python $(VENV)/bin/python -e '.[dev]')
	@echo "venv ready at $(VENV)"

key: ## Verify the .alphagenome_key file exists and looks plausible.
	@if [ ! -f $(KEY_FILE) ]; then echo "FAIL: $(KEY_FILE) missing"; exit 1; fi
	@LEN=$$(wc -c < $(KEY_FILE)); \
	if [ "$$LEN" -lt 30 ]; then echo "FAIL: $(KEY_FILE) is too short ($$LEN bytes)"; exit 1; fi
	@echo "OK: $(KEY_FILE) present ($$LEN bytes)"

check: ## Verify all expected outputs/*.csv exist (sanity check, no API calls).
	@missing=0; \
	for f in $(EXPECTED_OUTPUTS); do \
		if [ ! -f $$f ]; then echo "  MISSING: $$f"; missing=1; \
		else echo "  OK:      $$f"; fi; \
	done; \
	if [ $$missing -ne 0 ]; then echo "FAIL: $$missing expected outputs missing"; exit 1; fi
	@echo "PASS: all $(words $(EXPECTED_OUTPUTS)) expected outputs present"

smoke: venv key ## Re-score 5 SCN1A variants; verify within tolerance of recorded values (5 API calls).
	$(KEY_LAUNCHER) scripts/smoke_test.py

exp001: venv key ## Re-run Exp 001 ISM smoke (1 pathogenic + 1 benign, ~768 mutations total).
	ISM_SMOKE=1 $(KEY_LAUNCHER) research_notebook/experiments/001_ism_scn1a/ism_experiment.py

exp011: venv key ## Re-run Exp 011 ISM on the 4 Tier-1 SCN1A candidates (small).
	$(KEY_LAUNCHER) research_notebook/experiments/011_ism_tier1/ism_dnase_tier1.py

exp013: venv key ## Re-run Exp 013 Tier-1 cryptic-splice scoring (4 variants × 3 scorers).
	$(KEY_LAUNCHER) research_notebook/experiments/013_tier1_actionable/score_cryptic_splice.py

reproduce: check exp001 exp011 exp013 smoke  ## Run everything in order; fail fast.
	@echo
	@echo "=========================================="
	@echo "reproduce: all targets PASSED"
	@echo "=========================================="
	@echo "Outputs verified:"
	@for f in $(EXPECTED_OUTPUTS) $(EXP001_OUTPUTS) $(EXP011_OUTPUTS) $(EXP013_OUTPUTS); do \
		test -f $$f && echo "  ✓ $$f"; \
	done

test: venv ## Run the unit-test suite (does not require API key).
	$(VENV)/bin/pytest tests/ -v

clean: ## Remove generated test outputs (does not touch existing outputs/*.csv).
	rm -rf .pytest_cache tests/__pycache__
	@echo "cleaned: .pytest_cache, tests/__pycache__"
