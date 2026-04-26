# Development

The canonical contributing guide — quality standards, C++ style rules, architecture overview, and extension patterns — is [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md). This page summarises the key commands and CI structure.

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"   # builds the C++ extension + installs dev tools
```

The `.[dev]` extras install: `basedpyright`, `ruff`, `lizard`, `mkdocs-material`.

---

## Quality checks

All checks below must pass cleanly before merging.

### Lint and format (Python)

```bash
ruff check cfanalysis            # PEP8 (E), Pyflakes (F), isort (I) rules
ruff format --check cfanalysis   # formatting check; drop --check to auto-fix
```

Line length is 88. Run `ruff format cfanalysis` to auto-format.

### Type checking

```bash
python -m basedpyright cfanalysis
```

Must produce **0 errors and 0 warnings** — catches type mismatches and broken contracts before runtime, particularly at the `cfextract`/`cfanalysis` boundary where Python receives C++ structs. The suppressions in `pyproject.toml` (`reportMissingTypeStubs`, `reportUnknownMemberType`, etc.) are intentional — sklearn ships no type stubs. Do not remove them.

### Cyclomatic complexity

```bash
lizard src/ cfanalysis/ --CCN 10 --warnings_only -i -1
```

The threshold is 10 per function. The `-i -1` flag makes the check **informational** — it prints violations but does not fail CI. High cyclomatic complexity correlates with defect density and makes functions hard to test; aim to stay under 10. If a function genuinely needs more branches, document why.

### C++ build and tests

Use `/tmp` (not `build/`) to avoid CMake generator conflicts with the pip-driven build:

```bash
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml \
  --allow-running-no-tests
```

Tests use [Catch2 v3](https://github.com/catchorg/Catch2), fetched at configure time via CMake `FetchContent`.

The C++ build enables `-Wall -Wextra -Wpedantic` and a broad additional warning set. Do not suppress compiler warnings — in C++ they often indicate genuine bugs (uninitialised values, sign comparison, undefined behaviour paths), not just style issues.

---

## CMake options

| Option | Default | Description |
|---|---|---|
| `MAKE_TEST` | `OFF` | Build the Catch2 test binary (`test-cfextract`) |
| `MAKE_PY` | `OFF` | Build the pybind11 Python extension |

Both default to `OFF` in CMakeLists.txt. `pip install -e .` overrides `MAKE_PY=ON` via `pyproject.toml`. The manual cmake invocation above builds tests only, without the Python extension.

---

## CI

GitHub Actions (`.github/workflows/ci.yml`) runs three jobs:

**`checks`** — all non-tag pushes, PRs, and manual dispatches: lint (`ruff`), type check (`basedpyright`), complexity (`lizard`), C++ build + Catch2 tests, JUnit results published to GitHub Checks via [dorny/test-reporter](https://github.com/dorny/test-reporter).

**`docker_build_push`** — tag pushes (or manual dispatch for build-only): builds and pushes `ghcr.io/blex-max/als-challenge:<tag>` and `:latest` to GHCR.

**`docs`** — pushes to `main` only, after `checks` passes: runs `mkdocs gh-deploy --force` to rebuild and push the docs site to the `gh-pages` branch.

