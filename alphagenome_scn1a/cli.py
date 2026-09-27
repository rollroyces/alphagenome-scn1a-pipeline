"""Console entry point: `alphagenome-scn1a <subcommand>`.

Kept intentionally minimal — the real work lives in the per-experiment
scripts under research_notebook/experiments/* and scripts/*. This wrapper
exists so `pip install -e .` produces a usable CLI binary and so the
Makefile / Dockerfile can dispatch to it without hard-coding paths.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def _run(script: str, *args: str) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    print(f"[alphagenome-scn1a] {' '.join(cmd)}")
    return subprocess.call(cmd, cwd=str(REPO_ROOT))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="alphagenome-scn1a",
        description="CLI shim for the SCN1A / AlphaGenome reproducibility harness.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("smoke", help="Run scripts/smoke_test.py").set_defaults(
        handler=lambda _a: _run("scripts/smoke_test.py")
    )
    sub.add_parser("check", help="Verify outputs/*.csv exist (no API call)").set_defaults(
        handler=lambda _a: subprocess.call(["make", "check"], cwd=str(REPO_ROOT))
    )

    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    sys.exit(main())
