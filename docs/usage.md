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

```bash
git clone https://github.com/blex-max/als-challenge.git
cd als-challenge

python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

`pip install` triggers scikit-build-core, which compiles the C++ extension (`cfextract`) and installs the `cfanalysis` CLI. If `htslib` is not found automatically via `pkg-config`, pass the paths explicitly:

```bash
pip install . \
  -C cmake.define.HTSLIB_INCLUDE_DIR=/path/to/include \
  -C cmake.define.HTSLIB_LIBRARY=/path/to/libhts.a
```

On macOS with a non-default compiler (e.g. Apple Clang):
```bash
CC=/usr/bin/clang CXX=/usr/bin/clang++ pip install -e .
```

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
cfanalysis --manifest data/samples.csv --output-dir results/
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--manifest` | *(required)* | Path to the CSV manifest |
| `--contig` | `chr21` | Contig name to extract features from (must match BAM header) |
| `--output-dir` | `.` | Directory for all output files; created if absent |
| `--top-motifs` | `20` | Number of top-frequency 4-mer end motifs to include as classifier features |

### Output files

| File | Description |
|---|---|
| `classification_report.txt` | Sklearn LOO-CV classification report (precision, recall, F1 per class) |
| `sample_summary.csv` | Per-sample scalar statistics (fl_mean, methylation_mean, etc.) |
| `end_motifs.png` | Grouped bar chart of top-N end-motif frequencies, ALS vs CTRL |
| `frag_lengths.png` | Mean ± SEM fragment length distribution (0–600 bp) per group |
| `methylation.png` | Mean CpG methylation rate per group with individual sample points |
| `start_positions.png` | Mean read start position distribution across genomic bins |
| `end_positions.png` | Mean read end position distribution across genomic bins |

---

## Docker

A pre-built image is available on GHCR:

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
  --manifest /data/samples.csv \
  --output-dir /results
```

Build locally:

```bash
docker build -t cfanalysis .
```

