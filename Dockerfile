# Dockerfile for the SCN1A / AlphaGenome reproducibility harness.
# Image bundles the Python 3.13 venv, scripts/, data/, and pyproject.toml so
# that `docker run …` can hit the live AlphaGenome API and re-verify a small
# subset of the pipeline. The full pipeline is NOT re-run inside the image
# (that would burn the user's API quota).
#
# Build:
#     docker build -t alphagenome-scn1a:latest .
#
# Run (mount the API key file and the outputs dir so the smoke test can
# compare against recorded scores):
#     docker run --rm \
#         -v "$(pwd)/.alphagenome_key:/app/.alphagenome_key:ro" \
#         -v "$(pwd)/outputs:/app/outputs:ro" \
#         alphagenome-scn1a:latest
#
# The default CMD is `make smoke`. Override with e.g.
#     docker run --rm … alphagenome-scn1a:latest make check

FROM python:3.13-slim AS base

# System deps: only what's strictly needed to build numpy/pandas wheels if
# they're not already available as binaries, plus git + curl for uv install.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git \
    && rm -rf /var/lib/apt/lists/*

# Install uv (modern pip / venv replacement; matches the dev workflow).
RUN pip install --no-cache-dir uv

WORKDIR /app

# Copy project metadata first so dependency install can be cached.
COPY pyproject.toml ./

# Install the package + runtime deps. Use --system because we already have
# python:3.13-slim's Python and don't want a nested venv.
RUN uv pip install --system -e ".[dev]"

# Copy source needed for the reproducibility harness.
COPY Makefile ./
COPY scripts/ ./scripts/
COPY research_notebook/ ./research_notebook/
COPY tests/ ./tests/

# Outputs / data are mounted at runtime; nothing is baked into the image.
# The .alphagenome_key is also mounted at runtime (kept out of the image).
VOLUME ["/app/.alphagenome_key", "/app/outputs", "/app/data"]

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Default command: run the smoke test. Override with `docker run … <target>`.
CMD ["make", "smoke"]
