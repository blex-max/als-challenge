# cfDNA Analysis Pipeline

Cell-free DNA feature extraction and ALS vs Control classification from bisulfite-sequenced BAM files.

## Architecture

The pipeline is split into two layers:

```
cfextract   (C++ / pybind11)
│  Reads BAM files via htslib. Extracts per-fragment features and
│  reduces them to compact histograms and counts before returning
│  to Python. Memory scales with the number of distinct histogram
│  bins, not with the number of reads — important for 50 GB BAMs.
│
└── cfanalysis   (pure Python)
       Converts RegionMetrics into per-sample feature vectors,
       runs classification, and produces plots.
```

The key extension point is `cfanalysis/features.py::metrics_to_features()` — adding a new field to `RegionMetrics` requires one entry there; the rest of the pipeline picks it up automatically. See the [docs site](https://blex-max.github.io/als-challenge/) for the full data model, scalability analysis, and design rationale.

---

## Prerequisites

- CMake ≥ 3.22
- A C++17 compiler (see [Installation](#installation) note on macOS)
- htslib ≥ 1.14 — on macOS: `brew install htslib`
- Python ≥ 3.12

---

## Installation

Clone the repo, create a virtual environment, then pip-install (which builds the C++ extension via scikit-build-core):

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On macOS, explicitly use the Xcode system clang to avoid conflicts with any Homebrew LLVM on your `PATH`:

```bash
CC=/usr/bin/clang CXX=/usr/bin/clang++ pip install -e .
```

If htslib cannot be found automatically via pkg-config, pass the paths directly:

```bash
CC=/usr/bin/clang CXX=/usr/bin/clang++ pip install -e . \
  -C cmake.define.HTSLIB_INCLUDE_DIR=/path/to/include \
  -C cmake.define.HTSLIB_LIBRARY=/path/to/libhts.a
```

After any change to C++ source or `pyproject.toml`, re-run the same install command to rebuild.

---

## Local checks

With the virtual environment active and dev dependencies installed (`pip install -e ".[dev]"`):

```bash
# Format and lint
ruff check cfanalysis
ruff format --check cfanalysis   # drop --check to auto-fix

# Type checking
python -m basedpyright cfanalysis

# Cyclomatic complexity (informational — does not fail)
lizard src/ cfanalysis/ --CCN 10 --warnings_only -i -1

# C++ unit tests — use a /tmp build dir to avoid generator conflicts with the pip build
cmake -B /tmp/cfextract-build -DMAKE_TEST=ON -DMAKE_PY=OFF
cmake --build /tmp/cfextract-build --parallel
/tmp/cfextract-build/test-cfextract \
  --reporter console \
  --reporter junit::out=test-results.xml
```

These are the same steps run by the CI (`checks` job in `.github/workflows/ci.yml`).

---

## CI

Two jobs in `.github/workflows/ci.yml`:

| Job | Trigger |
|---|---|
| `checks` | All non-tag pushes and pull requests |
| `docker_build_push` | Tag push — builds and pushes Docker image to GHCR |

Both jobs can also be triggered manually via `workflow_dispatch`. To trigger and watch a run from the terminal (requires the [GitHub CLI](https://cli.github.com/)):

```bash
# Trigger — dispatches from whichever branch/ref you specify
gh workflow run CI --ref main

# Watch the latest run live
gh run watch

# Or pick a specific run
gh run list
gh run watch <run-id>
```

`docker_build_push` also runs on manual dispatch (build only — no push to GHCR).

---

## Usage

### Prepare a sample manifest

The pipeline takes a CSV with three columns:

```
sample_id,bam_path,label
SRR13404367,../data/bam/SRR13404367.bam,als
SRR13404371,../data/bam/SRR13404371.bam,ctrl
```

`bam_path` can be absolute or relative to the manifest file's directory. Each BAM must have a `.bai` index alongside it.

A manifest for the 22-sample ALS/CTRL cohort (chr21-only, downsampled) is at `data/samples.csv`.

### Run

```bash
cfanalysis --manifest data/samples.csv --output-dir results/
```

Options:

| Flag | Default | Description |
|---|---|---|
| `--manifest` | *(required)* | Path to sample manifest CSV |
| `--contig` | `chr21` | Contig name to process |
| `--output-dir` | `.` | Directory for output files |
| `--top-motifs` | `20` | Top-N end motifs to use as classification features |

### Output

| File | Description |
|---|---|
| `classification_report.txt` | LOO-CV precision / recall / F1 per class |
| `sample_summary.csv` | Per-sample scalar statistics (FL mean/std/short fraction, methylation mean/var/high fraction) |
| `end_motifs.png` | Grouped bar chart of top end-motif frequencies, ALS vs CTRL |
| `frag_lengths.png` | Mean ± SEM fragment length distribution per group (0–600 bp) |
| `methylation.png` | Bar chart of mean CpG methylation rate per group with individual sample points |
| `start_positions.png` | Fragment start position distribution per group (100 kbp bins) |
| `end_positions.png` | Fragment end position distribution per group (100 kbp bins) |

---

## Docker

Pull the pre-built image from GHCR:

```bash
docker pull ghcr.io/blex-max/als-challenge:latest
```

Or build locally:

```bash
docker build -t cfanalysis .
```

Then run, mounting your data (replace `<image>` with `ghcr.io/blex-max/als-challenge:latest` or `cfanalysis`):

```bash
docker run --rm \
  -v /path/to/bam:/data/bam \
  -v /path/to/manifest.csv:/data/samples.csv \
  -v /path/to/results:/results \
  <image> --manifest /data/samples.csv --output-dir /results
```

---

## Data

The 22-sample cohort is derived from [Li et al. (2021)](https://www.nature.com/articles/s41467-021-24843-4) (PRJNA691320): 12 ALS and 10 Control plasma cfDNA samples, bisulfite-sequenced on Illumina NovaSeq 6000 and aligned with Bismark. The working dataset is downsampled to 10 M reads per sample and restricted to chr21.

Full-genome BAMs are ~50 GB per sample; the pipeline is designed to handle this without loading all reads into memory (see [Scalability notes](#scalability-notes)).

---

## Current results (chr21)

LOO-CV with L2 logistic regression (C tuned per fold via inner 5-fold CV) over top-20 end motifs + fragment length statistics + CpG methylation statistics gives **0.64 accuracy / 0.62 macro F1** on the 22-sample cohort (12 ALS, 10 CTRL).

Fragment length features — particularly `fl_ratio_mono_di` (mono-nucleosomal / di-nucleosomal count ratio, effect size ~0.89) — carry most of the discriminative signal. CpG methylation features contribute minimally on chr21: the median CpG site is 100% methylated in every sample, and epipolymorphism (entropy) is identical across groups. Removing methylation features raises accuracy to 0.68 / macro F1 0.66 by reducing noise dimensions relative to n=22.

chr21 represents ~1.5% of the autosomal genome and contains no known ALS-associated differentially methylated regions. Performance is expected to increase substantially on full-genome BAMs.

