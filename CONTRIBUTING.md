# Contributing

## Development setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"   # builds the C++ extension via scikit-build-core + CMake
```

Note that pip installing should autogenerate typing stubs for the CPP core. These
should not usually need to be edited by hand.

The `.[dev]` extras install: `basedpyright`, `ruff`, `lizard`, `mkdocs-material`.

---

## Quality checks

All checks must pass before merging. CI runs them automatically on every push and PR; run locally before pushing.

### Lint and format (Python)

```bash
ruff check cfclassify            # PEP8 (E), Pyflakes (F), isort (I)
ruff format --check cfclassify   # drop --check to auto-fix
```

Line length is 88.

### Type checking

```bash
python -m basedpyright cfclassify
```

Must produce **0 errors and 0 warnings**. The suppressions in `pyproject.toml` (`reportMissingTypeStubs`, `reportUnknownMemberType`, etc.) are intentional — sklearn ships no type stubs. Do not remove them.

### Cyclomatic complexity

```bash
lizard src/ cfclassify/ --CCN 10 --warnings_only -i -1
```

Threshold is 10 per function. The `-i -1` flag is informational — violations are printed but CI does not fail. New code should conform; document the reason if a function genuinely needs more branches.

### C++ build and tests

Use `/tmp` to avoid CMake generator conflicts with the pip-driven build in `build/`:

```bash
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml
```

Tests use [Catch2 v3](https://github.com/catchorg/Catch2), fetched at configure time via CMake `FetchContent`. The build enables `-Wall -Wextra -Wpedantic` and a broad additional warning set — do not suppress compiler warnings to make code compile.

---

## C++ style

- Do not align expressions across multiple lines. The only exception is multiline inline comments.
- Function calls with arguments should have a space before the opening parenthesis.
- Function signatures: separate lines for return type, name, and argument list; opening brace on its own line.
- Non-trivial functions and classes should have a concise explanatory comment.

---

## Architecture

Two layers connected by the `RegionMetrics` boundary type:

```
cfextract   (C++ / pybind11)
│  src/core/access.{cpp,hpp}    — RAII BAM/index handle (AlnFile)
│  src/core/features.{cpp,hpp}  — per-read end-motif extraction
│  src/core/stats.{cpp,hpp}     — histogram → scalar summary stats
│  src/core/extract.{cpp,hpp}   — top-level extraction loop → RegionMetrics
│  src/bindings/python/         — pybind11 glue; exposes cfextract.extract_features()
│
│  ── RegionMetrics boundary ──
│
cfclassify   (pure Python)
   cfclassify/types.py      — TypedDict schema (Features, Sample)
   cfclassify/features.py   — metrics_to_features() — primary extension point
   cfclassify/classify.py   — LOO-CV with inner GridSearchCV
   cfclassify/plots.py      — matplotlib figures
   cfclassify/__main__.py   — CLI entry point
```

**Primary extension point**: adding a new feature to the pipeline requires changes to exactly two places:

1. Add a field to `RegionMetrics` in `src/core/extract.hpp` and populate it in `src/core/extract.cpp`.
2. Add the conversion logic to `cfclassify/features.py::metrics_to_features()`.

`build_feature_matrix()` in `classify.py` picks up any flat numeric key automatically — no further changes needed.

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

**`checks`** — all non-tag pushes, PRs, and manual dispatches: lint, type check, complexity, C++ build + Catch2 tests; JUnit XML published to GitHub Checks via `dorny/test-reporter`.

**`docker_build_push`** — tag pushes (or manual dispatch for build-only): builds and pushes `ghcr.io/blex-max/als-challenge:<tag>` and `:latest` to GHCR.

**`docs`** — pushes to `main` only, after `checks` passes: `mkdocs gh-deploy --force` rebuilds and deploys the docs site.

Branch model: feature branches → `main` directly.
