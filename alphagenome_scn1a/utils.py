"""API key handling, output formatting, and small helpers for the CLI.

Designed so every command goes through ``resolve_api_key`` and never prints
the actual key value — even partial keys can leak project quota.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path


# Where we look for the API key, in order:
#   1. $ALPHAGENOME_API_KEY (env var — preferred for CI / ephemeral shells)
#   2. ./.alphagenome_key    (in the current working directory)
#   3. ~/.alphagenome_key    (in the user's home directory)
#
# Anything else is rejected. We deliberately do NOT look in shell history,
# keyring, or other locations — keep it boring and reproducible.
_API_KEY_PATHS = (
    Path(".alphagenome_key"),
    Path.home() / ".alphagenome_key",
)

# AlphaGenome keys are ~40 chars and start with "AIzaSy". This is a heuristic
# (the format may evolve) but it catches the common "empty file" and "wrong
# file" mistakes at import time without round-tripping through the API.
_KEY_PATTERN = re.compile(r"^AIzaSy[A-Za-z0-9_\-]{30,}$")


@dataclass(frozen=True)
class ApiKey:
    """A validated API key. The value is in ``value``; never print it."""

    value: str
    source: str  # "env" | ".alphagenome_key" | "~/.alphagenome_key"

    def masked(self) -> str:
        """Return a masked form like ``AIzaSy****xyz`` — safe to print."""
        if len(self.value) <= 8:
            return "****"
        return f"{self.value[:6]}****{self.value[-4:]}"

    def is_plausible(self) -> bool:
        """True if the key passes our format heuristic."""
        return bool(_KEY_PATTERN.match(self.value))


def resolve_api_key(require: bool = True) -> ApiKey | None:
    """Find the API key, validate its format, and return it.

    Args:
        require: when True, raise ``MissingApiKeyError`` if nothing is found.
                 When False, return None instead — useful for commands that
                 don't actually need the key (``info``, ``tier1``).

    Raises:
        MissingApiKeyError: the key wasn't found anywhere (only if require=True).
        InvalidApiKeyError: the key was found but doesn't look like a real one.
    """
    raw = os.environ.get("ALPHAGENOME_API_KEY")
    if raw:
        return _validate(raw, source="env")

    for path in _API_KEY_PATHS:
        if path.exists():
            text = path.read_text().strip()
            if text:
                return _validate(text, source=str(path))

    if require:
        raise MissingApiKeyError()
    return None


def api_key_status_lines() -> list[str]:
    """Return formatted lines for the ``info`` command.

    Three possible states:
      - resolved and format-valid → one masked-key line
      - resolved but doesn't match format → one warning line
      - not found at all → one "missing" line
    """
    try:
        key = resolve_api_key(require=False)
    except InvalidApiKeyError as e:
        return [f"warning — {e}"]

    if key is None:
        return [
            "status:       NOT FOUND",
            "expected:     $ALPHAGENOME_API_KEY or ./.alphagenome_key or ~/.alphagenome_key",
        ]
    valid = "yes" if key.is_plausible() else "no (wrong format)"
    return [
        f"status:       resolved ({key.source})",
        f"format check: {valid}",
        f"masked:       {key.masked()}",
    ]


def _validate(value: str, source: str) -> ApiKey:
    """Internal — check the value's shape before returning it."""
    if not value or not value.strip():
        raise InvalidApiKeyError(
            "API key file is empty.",
            source=source,
        )
    if len(value) < 30:
        raise InvalidApiKeyError(
            f"API key is only {len(value)} chars — expected ~40.",
            source=source,
        )
    return ApiKey(value=value, source=source)


class MissingApiKeyError(RuntimeError):
    """Raised when no API key can be located anywhere."""

    def __init__(self) -> None:
        super().__init__(
            "No AlphaGenome API key found.\n"
            "\n"
            "  How to fix:\n"
            "    1. Get a free key at https://alphagenome.google/api\n"
            "    2. Either:\n"
            "         export ALPHAGENOME_API_KEY=AIzaSy...\n"
            "       or:\n"
            "         echo 'AIzaSy...' > .alphagenome_key && chmod 600 .alphagenome_key\n"
            "\n"
            "  Commands that DO need a key: score-gene, score-vcf, reproduce\n"
            "  Commands that do NOT:        info, tier1\n"
        )


class InvalidApiKeyError(RuntimeError):
    """Raised when the key was found but doesn't look right."""

    def __init__(self, message: str, source: str) -> None:
        super().__init__(
            f"{message} (source: {source}). Get a fresh key at "
            f"https://alphagenome.google/api"
        )
        self.source = source


def humanize_path(path: str | Path) -> str:
    """Pretty-print a path relative to CWD when possible."""
    p = Path(path).resolve()
    try:
        return str(p.relative_to(Path.cwd()))
    except ValueError:
        return str(p)


def banner(title: str, char: str = "=") -> str:
    """Return a 70-char centered banner. Used to make CLI output scannable."""
    width = 70
    if len(title) >= width - 4:
        return f" {title} "
    pad = (width - len(title) - 2) // 2
    return f"{char * pad} {title} {char * pad}"


def fmt_count(n: int, singular: str, plural: str | None = None) -> str:
    """``1 variant`` vs ``3 variants``."""
    if n == 1:
        return f"{n} {singular}"
    return f"{n} {plural or (singular + 's')}"