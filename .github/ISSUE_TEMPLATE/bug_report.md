---
name: Bug report
about: Report something that isn't working (broken output, failing test, etc.)
title: "[bug] "
labels: ["bug"]
assignees: []
---

### What broke

A clear one-paragraph description of what went wrong.

### How to reproduce

Steps to reproduce the behaviour:

```bash
# exact commands you ran
```

### What I expected

What you expected to happen.

### What actually happened

```
# paste the actual error / output / log here
```

### Environment

- OS / runner:
- Python version (`python3 --version`):
- `alphagenome` SDK version (`uv pip show alphagenome`):
- Branch / commit (`git rev-parse HEAD`):
- API key present? (`[ -f .alphagenome_key ] && echo yes || echo no`):

### Relevant files

- Output file (e.g. `outputs/vus_rescored.csv`):
- Script that failed (e.g. `scripts/rescore_vus.py`):

### Additional context

Anything else that might help — screenshots, links, related issues/PRs.