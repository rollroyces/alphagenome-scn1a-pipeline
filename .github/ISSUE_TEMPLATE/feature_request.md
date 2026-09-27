---
name: Feature request
about: Suggest a new analysis, output, or workflow improvement
title: "[feature] "
labels: ["enhancement"]
assignees: []
---

### Problem / motivation

What limitation or gap in the current pipeline prompted the idea? Quote a
specific shortcoming if you can — e.g. "the VUS table has no DNase fold-change
column" or "we have no way to compare against a second model".

### Proposed solution

Describe the change you'd like, ideally as a small, testable delta:

- New output file (path + format):
- New / changed script (entry point + args):
- New Makefile target (if any):

### Alternatives considered

Any other approaches you weighed, and why you prefer the one above.

### Acceptance criteria

How would we know this is done? E.g.

- [ ] `make exp013` produces `outputs/<new_file>.csv`
- [ ] new pytest under `tests/` covers the changed function
- [ ] README's Quickstart mentions the new entry point

### Out of scope

Anything you explicitly want to defer to a follow-up issue.