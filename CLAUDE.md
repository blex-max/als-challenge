# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Development setup
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"   # builds the C++ extension via scikit-build-core + CMake
```

### Lint and format (Python)
```bash
ruff check cfanalysis            # lint (E, F, I rules)
ruff format --check cfanalysis   # format check; drop --check to auto-fix
python -m basedpyright cfanalysis
```

### Complexity (C++ and Python)
```bash
lizard src/ cfanalysis/ --CCN 10 --warnings_only -i -1
```

### C++ build and tests
```bash
# Use /tmp to avoid CMake generator conflicts with the pip-driven build in build/
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml \
  --allow-running-no-tests
```

## Architecture

The pipeline has two layers connected by the `RegionMetrics` boundary type:

```
cfextract   (C++ / pybind11)
│  Reads BAM files via htslib. Extracts per-fragment features and reduces
│  them to compact histograms before returning to Python. Memory is O(bins),
│  not O(reads) — important for 50 GB BAMs.
│
└── cfanalysis   (pure Python)
       Converts RegionMetrics into per-sample feature vectors,
       runs LOO-CV classification, and produces plots.
```

Key source locations:

| Path | Role |
|---|---|
| `src/core/` | C++ extraction logic (access, extract, features) |
| `src/bindings/python/bindings.cpp` | pybind11 glue — exposes `cfextract.extract_features()` |
| `cfanalysis/features.py` | `metrics_to_features()` — the primary extension point |
| `cfanalysis/classify.py` | LOO-CV SVM classification |
| `cfanalysis/plots.py` | End-motif bar chart |
| `tests/cpp/test_extract.cpp` | Catch2 v3 C++ tests (fetched via FetchContent) |

**Primary extension point:** adding a new field to `RegionMetrics` (C++ struct) requires a corresponding entry in `cfanalysis/features.py::metrics_to_features()` — the rest of the pipeline picks it up automatically.

There are no Python-level tests yet; the test infrastructure is C++-only.

## Quality standards

- **Cyclomatic complexity:** max 10 per function (lizard, `--CCN 10`). The check is informational (`-i -1`) and does not fail the build.
- **Type checking:** `basedpyright cfanalysis` must produce 0 errors and 0 warnings. The sklearn-related suppressions in `pyproject.toml` (`reportMissingTypeStubs`, `reportUnknownMemberType`, etc.) are intentional — sklearn ships no type stubs. Do not remove them.
- **Format and lint:** ruff with E, F, I rules and line length 88. Both `ruff check` and `ruff format --check` must pass cleanly.
- **C++ warnings:** the CMake build enables `-Wall -Wextra -Wpedantic` and a broad set of additional flags. Do not suppress compiler warnings to make code compile.

## CI

GitHub Actions (`.github/workflows/ci.yml`), two jobs:

- **`checks`** — all non-tag events (pushes, PRs, manual dispatch): lint, type check, complexity, C++ build + tests; JUnit test results uploaded as artifact
- **`docker_build_push`** — tag pushes: builds and pushes Docker image to GHCR (`ghcr.io/<repo>:<tag>` and `:latest`); also runs on manual dispatch (build only, no push)

Test results are published to GitHub Checks via `dorny/test-reporter` (JUnit XML from Catch2).

Branch model: feature branches → `main` directly (no develop branch).
