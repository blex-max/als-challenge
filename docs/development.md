# Development

The canonical contributing guide is [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md). This page explains the rationale for those standards.

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
ruff check cfclassify            # PEP8 (E), Pyflakes (F), isort (I) rules
ruff format --check cfclassify   # formatting check; drop --check to auto-fix
```

Line length is 88. Run `ruff format cfclassify` to auto-format.

A consistent style eliminates formatting noise from diffs and makes code review faster. Enforcing it in CI means contributors never have to debate formatting choices. Ruff is used over flake8+black because it covers linting, formatting, and import ordering in a single fast pass.

### Type checking

```bash
python -m basedpyright cfclassify
```

Must produce **0 errors and 0 warnings** — catches type mismatches and broken contracts before runtime, particularly at the `cfextract`/`cfclassify` boundary where Python receives C++ structs. The suppressions in `pyproject.toml` (`reportMissingTypeStubs`, `reportUnknownMemberType`, etc.) are intentional since unfortunately sklearn ships no type stubs.
 basedpyright is used over mypy/pyright because it enforces stricter defaults.

### Cyclomatic complexity

```bash
lizard src/ cfclassify/ --CCN 10 --warnings_only -i -1
```

lizard reports a cross-language assessment of code complexity. High cyclomatic complexity (CCN) correlates with defect density and makes functions harder to reason about both for humans and AI coding tools. The informational threshold (-i) guides design without penalising genuinely complex logic: the goal is to flag cases worth reconsidering, not to enforce a hard rule that forces artificial splitting.

### C++ build and tests

```bash
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml
```

Tests use [Catch2 v3](https://github.com/catchorg/Catch2), fetched at configure time via CMake `FetchContent`.

Testing at the unit level documents the intended behaviour of each statistical computation independently of the extraction loop, making regressions and bugs easy to locate.

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

