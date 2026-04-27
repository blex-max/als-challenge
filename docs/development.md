# Development

The canonical contributing guide — quality standards, C++ style rules, architecture overview, and extension patterns — is [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md). This page explains the *rationale* behind those standards. Commands are included as a reference, but the emphasis is on *why* each standard exists — because a future contributor (or an interviewer) should understand the intent, not just the mechanics.

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

**Rationale**: a consistent style eliminates formatting noise from diffs and makes code review faster. Enforcing it in CI means contributors never debate formatting choices — the tool is authoritative. Ruff is used over flake8+black because it covers linting, formatting, and import ordering in a single fast pass.

### Type checking

```bash
python -m basedpyright cfclassify
```

Must produce **0 errors and 0 warnings** — catches type mismatches and broken contracts before runtime, particularly at the `cfextract`/`cfclassify` boundary where Python receives C++ structs. The suppressions in `pyproject.toml` (`reportMissingTypeStubs`, `reportUnknownMemberType`, etc.) are intentional — sklearn ships no type stubs. Do not remove them.

**Rationale**: the `cfextract`/`cfclassify` boundary is where pybind11 objects (with no Python stubs by default) meet scikit-learn types. Without strict type checking, `AttributeError` and shape mismatches at this boundary are only caught at runtime on real BAM data. The 0-warnings policy keeps the type graph complete — no implicit `Any` holes that could mask a broken contract after an API change. basedpyright is used over mypy because it enforces stricter defaults and does not silently infer `Unknown` as `Any`.

### Cyclomatic complexity

```bash
lizard src/ cfclassify/ --CCN 10 --warnings_only -i -1
```

The threshold is 10 per function. The `-i -1` flag makes the check **informational** — it prints violations but does not fail CI. High cyclomatic complexity correlates with defect density and makes functions hard to test; aim to stay under 10. If a function genuinely needs more branches, document why.

**Rationale**: high cyclomatic complexity (CCN) correlates with defect density and makes functions harder to reason about — both for humans and AI coding tools. The informational threshold (not a CI blocker) guides design without penalising genuinely complex logic: the goal is to flag cases worth reconsidering, not to enforce a hard rule that forces artificial splitting.

### C++ build and tests

Use `/tmp` (not `build/`) to avoid CMake generator conflicts with the pip-driven build:

```bash
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml
```

Tests use [Catch2 v3](https://github.com/catchorg/Catch2), fetched at configure time via CMake `FetchContent`.

The C++ build enables `-Wall -Wextra -Wpedantic` and a broad additional warning set. Do not suppress compiler warnings — in C++ they often indicate genuine bugs (uninitialised values, sign comparison, undefined behaviour paths), not just style issues.

**Rationale**: C++ compiler warnings frequently indicate real defects — uninitialised values, sign conversion mistakes, shadowed variables — not stylistic preferences. Suppressing them with `#pragma GCC diagnostic ignore` or `-Wno-*` trades a compile-time diagnostic for a potential runtime bug. The `-Wall -Wextra -Wpedantic` set is the standard minimum for production C++ and is a prerequisite for sanitisers (ASAN, UBSAN) to give meaningful output.

**Rationale (Catch2 tests)**: the unit tests cover edge cases (empty histograms, fully methylated cohorts, zero-fragment inputs) that are impractical to exercise through integration testing — constructing pathological BAMs just to hit these branches would be fragile and slow. Testing at the unit level also documents the intended behaviour of each statistical computation independently of the extraction loop, making regressions immediately locatable.

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

