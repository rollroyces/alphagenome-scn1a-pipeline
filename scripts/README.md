# Scripts

Reusable scripts for AlphaGenome exploration. Run them with the venv active
and `ALPHAGENOME_API_KEY` set.

## Smoke tests (run first)

| Script | What it verifies | Runtime | Quota cost |
|---|---|---|---|
| `smoke_test.py` | API key works, model endpoint reachable, trivial prediction round-trips | ~30s | 1 prediction |
| `atlas_smoke_test.py` | Atlas access works, returns pre-computed scores | ~1 min | 1 interval query |

Run smoke tests *before* opening the notebooks so failures are easy to debug.

## Examples

| Script | Equivalent notebook | Purpose |
|---|---|---|
| `quickstart_script.py` | `quick_start.ipynb` (first 5 sections) | Predict / score / variant effect in one shot |
| `ism_example.py` | `quick_start.ipynb` (ISM section) | In-silico mutagenesis walkthrough |

These exist because Jupyter adds friction when you just want to "does this
work?" — pure scripts give faster feedback and can be diff'd cleanly.

## Conventions

- All scripts print a `=== Step N ===` header at the start of each phase
- All scripts write figures to `figures/`
- All scripts exit non-zero on error so they work in CI / cron
- All scripts accept `ALPHAGENOME_API_KEY` from env, never as a CLI arg
