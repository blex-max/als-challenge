# Usage

## Prerequisites

| Dependency | Minimum version | Notes |
|---|---|---|
| CMake | 3.22 | Required to build the C++ extension |
| C++17 compiler | — | clang++ or g++ |
| htslib | 1.14 | BAM/CRAM I/O; install via package manager or from source |
| Python | 3.12 | |

On macOS (Homebrew):
```bash
brew install cmake htslib
```

On Ubuntu/Debian:
```bash
sudo apt-get install cmake build-essential pkg-config libhts-dev
```

---

## Installation

### Option A — Docker (recommended for reproducibility)

A pre-built image is available on GHCR and bundles all C++ dependencies. No local build is needed.

```bash
docker pull ghcr.io/blex-max/als-challenge:latest
```

Run the pipeline with volume mounts:

```bash
docker run --rm \
  -v /host/path/to/bams:/data/bams \
  -v /host/path/to/samples.csv:/data/samples.csv \
  -v /host/path/to/results:/results \
  ghcr.io/blex-max/als-challenge:latest \
  train --manifest /data/samples.csv --out-dir /results
```

Build locally:

```bash
docker build -t cfclassify .
```

### Option B — From source

Before running, prepare a manifest CSV describing your samples — see [Manifest file](#manifest-file) below.

```bash
git clone https://github.com/blex-max/als-challenge.git
cd als-challenge

python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

`pip install` triggers scikit-build-core, which compiles the C++ extension (`cfextract`) and installs the `cfclassify` CLI. If `htslib` is not found automatically via `pkg-config`, pass the paths explicitly:

```bash
pip install . \
  -C cmake.define.HTSLIB_INCLUDE_DIR=/path/to/include \
  -C cmake.define.HTSLIB_LIBRARY=/path/to/libhts.a
```

If Homebrew's Clang takes precedence over the system compiler (common on macOS), specify the compiler explicitly:

```bash
CC=/usr/bin/clang CXX=/usr/bin/clang++ pip install -e .
```

This is typically needed on Apple Silicon Macs where Xcode CLT and Homebrew Clang coexist.

---

## Manifest file

The pipeline takes a CSV manifest describing the samples:

```csv
sample_id,bam_path,label
ALS_001,bams/ALS_001.bam,als
ALS_002,bams/ALS_002.bam,als
CTRL_001,bams/CTRL_001.bam,ctrl
```

- `sample_id` — unique identifier, used in plot labels and the summary CSV
- `bam_path` — path to the BAM file; relative paths are resolved from the directory containing the manifest
- `label` — group label; the classifier treats these as class names (any two distinct strings are valid)

Each BAM file must have an accompanying BAI index at `<bam>.bai` or `<bam>.csi`.

---

## Running

```bash
cfclassify <subcommand> [options]
```

Three subcommands are available.

### `train` — Train on a labelled cohort

```bash
cfclassify train --manifest data/samples.csv --out-dir results/
```

Runs LOO-CV and writes the evaluation report and plots. Also fits a final model on all samples and saves it to `<out-dir>/model.pkl` alongside a feature cache for use by `update`.

| Flag | Default | Description |
|---|---|---|
| `--manifest` | *(required)* | Path to the CSV manifest |
| `--contig` | `chr21` | Contig name to extract features from (must match BAM header) |
| `--out-dir` | `.` | Directory for all output files; created if absent |
| `--motif-length` | `4` | k-mer length for end-motif features; 4^k features used (default 4 → 256) |

Output files:

| File | Description |
|---|---|
| `classification_report.txt` | Sklearn LOO-CV classification report (precision, recall, F1 per class) |
| `sample_summary.csv` | Per-sample scalar statistics (fl_mean, methylation_mean, etc.) |
| `model.pkl` | Trained model bundle (joblib) |
| `model.pkl.features.json` | Feature cache for `update` |
| `end_motifs.png` | Grouped bar chart of top-10 end-motif frequencies, ALS vs CTRL |
| `frag_lengths.png` | Mean ± SEM fragment length distribution (0–600 bp) per group |
| `methylation.png` | Mean CpG methylation rate per group with individual sample points |
| `start_positions.png` | Mean read start position distribution across genomic bins |
| `end_positions.png` | Mean read end position distribution across genomic bins |

### `predict` — Classify a single unlabelled BAM

```bash
cfclassify predict \
  --bam patient_001.bam \
  --model-path models/als_chr21.pkl \
  --sample-id patient_001
```

| Flag | Default | Description |
|---|---|---|
| `--bam` | *(required)* | Path to the BAM file (must be indexed) |
| `--model-path` | *(required)* | Path to a saved model bundle |
| `--sample-id` | `unknown` | Label for progress output |

Prints the predicted class label and probability to stdout. The contig and k-mer length are read from the saved bundle — no separate flags needed.

### `update` — Add a new labelled sample and retrain

```bash
cfclassify update \
  --bam new_ctrl.bam \
  --label ctrl \
  --sample-id CTRL_023 \
  --model-path models/als_chr21.pkl
```

| Flag | Default | Description |
|---|---|---|
| `--bam` | *(required)* | BAM file for the new sample |
| `--label` | *(required)* | Class label for the new sample |
| `--sample-id` | *(required)* | Unique identifier |
| `--model-path` | *(required)* | Path to an existing model bundle |

Extracts features using the stored contig and k, appends the sample to the feature cache, and retrains the model from scratch on all cached samples. Updates both `<model-path>` and `<model-path>.features.json` in place.

Requires at least 2 samples per class in the cache. Run `train` separately if updated LOO-CV metrics are needed.


