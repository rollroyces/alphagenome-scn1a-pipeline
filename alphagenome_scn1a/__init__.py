"""Public package surface for ``alphagenome-scn1a``.

The package exists primarily so ``pip install -e .`` produces a usable
``alphagenome-scn1a`` console script (declared in ``pyproject.toml``).

The real pipeline logic lives in :mod:`alphagenome_scn1a.cli` (Typer
commands) and the per-experiment scripts under ``scripts/``. The
modules here (``extract``, ``score``, ``pipeline``, ``utils``) are
importable helpers used by the CLI; nothing here is private to the CLI.
"""

from __future__ import annotations

__version__ = "0.3.0"

__all__ = [
    "__version__",
    "extract",
    "score",
    "pipeline",
    "utils",
]


# Re-export the submodules so callers can do
# ``from alphagenome_scn1a import extract, score, pipeline``.
from . import extract, score, pipeline, utils  # noqa: E402, F401