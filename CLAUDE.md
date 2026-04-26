# CLAUDE.md

Guidance for Claude Code working in this repository. See [CONTRIBUTING.md](CONTRIBUTING.md) for the full development guide: commands, quality standards, C++ style rules, architecture overview, and CI.

## Key reminders

- Run `lizard` when verifying any work — CCN ≤ 10 per function.
- `basedpyright cfanalysis` must produce 0 errors/warnings. The sklearn suppressions in `pyproject.toml` are intentional; do not remove them.
- C++ test builds: use `/tmp/cfextract-build`, not `build/`, to avoid generator conflicts with the pip-driven build.
- Do not suppress compiler warnings in C++ to make code compile.
