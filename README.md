# ALS cfDNA Analysis

**FULL GUIDE AND DOCUMENTATION AVAILABLE [HERE](https://blex-max.github.io/als-challenge/)**

cfDNA feature-extraction and exploratory classification pipeline for distinguishing ALS and control samples using bisulfite-sequenced plasma BAMs.

22-sample chr21 cohort · 266 features · **0.64 accuracy / 0.62 macro F1** under LOO-CV.

---

## Quick start

Create a manifest CSV listing your samples:

```csv
sample_id,bam_path,label
ALS_001,bams/ALS_001.bam,als
CTRL_001,bams/CTRL_001.bam,ctrl
```

**Option A — Docker (no local build required)**

```bash
docker pull ghcr.io/blex-max/als-challenge:latest

# Train: runs LOO-CV evaluation and saves a deployable model to /results/model.pkl
docker run --rm \
  -v /path/to/bams:/data/bams \
  -v /path/to/samples.csv:/data/samples.csv \
  -v /path/to/results:/results \
  ghcr.io/blex-max/als-challenge:latest \
  train --manifest /data/samples.csv --out-dir /results

# Predict: classify a new unlabelled sample against the saved model
docker run --rm \
  -v /path/to/bams:/data/bams \
  -v /path/to/results:/results \
  ghcr.io/blex-max/als-challenge:latest \
  predict --bam /data/bams/new_patient.bam \
          --model-path /results/model.pkl \
          --sample-id new_patient
```

**Option B — pip install (from source)**

```bash
git clone https://github.com/blex-max/als-challenge.git && cd als-challenge
python3 -m venv .venv && source .venv/bin/activate
pip install -e .

cfclassify train --manifest data/samples.csv --out-dir results/
cfclassify predict --bam new_patient.bam --model-path results/model.pkl --sample-id new_patient
```

---

## Features

| Feature class | Key signal |
|---|---|
| **End-motif frequencies** | Differential appearance of certain k-mers at fragment ends |
| **Fragment length** | Length distribution of cfDNA fragments |
| **CpG methylation** | Differential methylation between samples |

---

## Documentation

- [Architecture](docs/architecture.md) — C++/Python boundary, memory and time complexity, extension patterns
- [Classification design](docs/classification.md) — feature engineering, LOO-CV, model persistence
- [Results](docs/results.md) — plots, per-sample data, and classification metrics
- [Usage / CLI reference](docs/usage.md) — prerequisites, installation, subcommand reference
- [Python API](docs/api/python.md) — `cfclassify` module reference
- [C++ API](docs/api/cpp.md) — `cfextract` library reference
- [Contributing](CONTRIBUTING.md) — development setup, quality checks, C++ style

---

## References

1. Caggiano C, Celona B, Garton F, et al. *Comprehensive cell type decomposition of circulating cell-free DNA with CelFiE*. Nature Communications. 2021;12:2717. https://doi.org/10.1038/s41467-021-22901-x

2. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. *Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin*. Cell. 2016;164(1-2):57–68. https://doi.org/10.1016/j.cell.2015.11.050

3. Ding SC, Lo YMD. *Cell-Free DNA Fragmentomics in Liquid Biopsy*. Diagnostics. 2022;12(4):978. https://doi.org/10.3390/diagnostics12040978

4. Moss J, Magenheim J, Neiman D, et al. *Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease*. Nature Communications. 2018;9:5068. https://doi.org/10.1038/s41467-018-07466-6
