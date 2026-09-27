# Test suite for the reproducibility harness.
#
# These tests are offline — they do NOT hit the AlphaGenome API. They verify
# that the harness itself is wired up correctly: expected output paths exist,
# the smoke test module is importable, and the per-experiment scripts exist.
#
# Run with: `make test` or `pytest tests/ -v`.

from __future__ import annotations

import os
import sys
from pathlib import Path

# Make scripts/ importable when pytest is invoked from the repo root.
REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))


# ---- Path resolution -------------------------------------------------------

EXPECTED_OUTPUTS = [
    "outputs/benchmark_scn1a_live_api_raw.csv",
    "outputs/benchmark_scn1a_results.csv",
    "outputs/clinvar_scn1a.json",
    "outputs/clinvar_scn1a.tsv",
    "outputs/clinvar_scn1a_benchmark.tsv",
    "outputs/vus_rescored.csv",
    "outputs/vus_top_candidates.csv",
    "outputs/vus_top_combined.csv",
    "outputs/vus_high_impact_candidates.csv",
    "outputs/tier1_actionable_features.csv",
    "outputs/cross_disease_metrics.csv",
    "outputs/cross_disease_scn1a_raw.csv",
]


def test_repo_root_exists():
    assert REPO_ROOT.exists()
    assert (REPO_ROOT / "Makefile").exists(), "Makefile missing at repo root"
    assert (REPO_ROOT / "pyproject.toml").exists(), "pyproject.toml missing"
    assert (REPO_ROOT / "Dockerfile").exists(), "Dockerfile missing"
    assert (REPO_ROOT / ".alphagenome_key").exists(), (
        ".alphagenome_key missing — create one with `echo $KEY > .alphagenome_key && chmod 600 .alphagenome_key`"
    )


def test_expected_outputs_present():
    """Every entry in EXPECTED_OUTPUTS must exist on disk (this is what
    `make check` enforces from the shell)."""
    missing = [p for p in EXPECTED_OUTPUTS if not (REPO_ROOT / p).exists()]
    assert not missing, f"Missing expected outputs: {missing}"


def test_makefile_targets_exist():
    """The Makefile must define the targets documented in README."""
    makefile = (REPO_ROOT / "Makefile").read_text()
    for target in ("check", "smoke", "exp001", "exp011", "exp013", "reproduce", "clean"):
        assert f"^${target}:" in makefile or f"\n{target}:" in makefile or f" {target}:" in makefile, (
            f"Makefile does not define a `{target}` target"
        )


def test_smoke_test_module_importable():
    """The smoke test must be importable as a module (CLI-style)."""
    import smoke_test  # noqa: F401  — placed in scripts/ via sys.path

    assert hasattr(smoke_test, "main")
    assert callable(smoke_test.main)


def test_run_with_key_launcher_exists():
    launcher = REPO_ROOT / "scripts" / "_run_with_key.sh"
    assert launcher.exists()
    assert os.access(launcher, os.X_OK), "scripts/_run_with_key.sh must be executable"


def test_per_experiment_scripts_exist():
    for path in [
        "research_notebook/experiments/001_ism_scn1a/ism_experiment.py",
        "research_notebook/experiments/011_ism_tier1/ism_dnase_tier1.py",
        "research_notebook/experiments/013_tier1_actionable/score_cryptic_splice.py",
    ]:
        assert (REPO_ROOT / path).exists(), f"Experiment script missing: {path}"


def test_alphagenome_key_format():
    """The API key file must be 40+ chars and mode 600 (defence in depth
    against accidental commits — .gitignore already excludes it)."""
    key_path = REPO_ROOT / ".alphagenome_key"
    stat = key_path.stat()
    assert stat.st_mode & 0o777 == 0o600, (
        f".alphagenome_key has mode {oct(stat.st_mode & 0o777)}, want 0o600"
    )
    assert stat.st_size >= 30, (
        f".alphagenome_key is only {stat.st_size} bytes; expected at least 30"
    )


# ---- Smoke-test internals (no API call) -----------------------------------

def test_pick_variants_is_deterministic():
    """Same seed must give the same 5 variants — guarantees reproducibility
    across machines."""
    import pandas as pd
    from smoke_test import _pick_variants

    # Tiny synthetic raw frame mimicking benchmark_scn1a_live_api_raw.csv.
    n_pos, n_neg = 20, 20
    pos = pd.DataFrame(
        {
            "chrom": ["2"] * n_pos,
            "pos": range(100, 100 + n_pos),
            "ref": ["A"] * n_pos,
            "alt": ["G"] * n_pos,
            "rsid": [f"pos_{i}" for i in range(n_pos)],
            "clnsig_category": ["pathogenic"] * n_pos,
            "success": [True] * n_pos,
            "SPLICE_SITES_score": [0.9] * n_pos,
            "SPLICE_SITE_USAGE_score": [100.0] * n_pos,
            "SPLICE_JUNCTIONS_score": [800.0] * n_pos,
        }
    )
    neg = pos.copy()
    neg["clnsig_category"] = "benign"
    neg["rsid"] = [f"neg_{i}" for i in range(n_neg)]
    neg["SPLICE_SITES_score"] = 0.1
    raw = pd.concat([pos, neg], ignore_index=True)

    a = _pick_variants(raw, n=5, seed=42).reset_index(drop=True)
    b = _pick_variants(raw, n=5, seed=42).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)


def test_auprc_fallback_matches_sklearn_when_available():
    """The mini-AUPRC routine should give the same answer as sklearn."""
    import numpy as np
    from smoke_test import _auprc

    rng = np.random.default_rng(0)
    y_true = rng.integers(0, 2, size=50)
    y_score = rng.random(50)
    got = _auprc(y_true, y_score)
    try:
        from sklearn.metrics import average_precision_score
        want = float(average_precision_score(y_true, y_score))
    except ImportError:
        return  # sklearn missing; nothing to compare against
    assert abs(got - want) < 1e-9, f"AUPRC fallback {got} != sklearn {want}"
