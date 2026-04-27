# ALS cfDNA Analysis

A cfDNA feature-extraction and classification pipeline for identifying putative ALS patients using bisulfite-sequenced plasma cell-free DNA.

---

## Background

Cell-free DNA (cfDNA) in blood plasma is shed by apoptotic and necrotic cells throughout the body. Because cfDNA carries epigenetic and fragmentation signatures from its tissue of origin, it can serve as a liquid biopsy for a range of conditions — including neurodegenerative diseases such as ALS.

This tool addresses the question: can we distinguish ALS patients from healthy controls using features extracted from bisulfite-sequenced plasma cfDNA BAM files?

---

## Pipeline at a glance

The C++ core (`cfextract`) handles BAM I/O via htslib and accumulates per-read statistics (end motifs, fragment lengths, and CpG methylation) into a compact `RegionMetrics` boundary object. The Python layer (`cfclassify`) converts that boundary object to a flat feature vector and trains or applies an L2 logistic regression classifier. On the chr21 cohort of 22 bisulfite-sequenced plasma BAMs (12 ALS, 10 CTRL), the pipeline extracts 266 features per sample and reaches 0.64 accuracy under LOO-CV.

---

## Engineering highlights

- **CI on every push** — a series of code analysis tools, and unit tests, run on every commit; failures surface in GitHub Checks, so regressions and bugs are caught before they reach a production environment.
- **Portable container** — pre-built image available via package registry (`docker pull ghcr.io/blex-max/als-challenge:latest`) eliminates end-user build issues; the full pipeline runs reproducibly on any machine with a single pull.
- **Language-agnostic extraction core** — core feature extraction and file I/O are handled in performant cpp. Python bindings are provided out of the box, but analysts can drive analysis from Python, R, Julia, or any language for which bindings can be made, without having to reimplement key functionality.
- **Memory scales with region size only, not read depth** — peak footprint is bounded by the number of unique CpG sites in the target region, not read count; a deeply sequenced chr21 BAM uses the same memory as a shallow one, keeping extraction tractable on standard hardware.
- **Fast extraction** — `cfextract.extract_features()` completes chr21 extraction in **0.10 ± 0.00 s** wall time, 5.9 ± 0.3 MB peak RSS across 22 chr21 BAMs (10 M reads per sample, Apple M-series, 3 repeats; reproduce with `python scripts/bench_extract.py --repeat 3`).
- **Deployable model** — the `train` subcommand fits and persists a final model bundle; `predict` applies it to new samples without retraining, enabling deployment beyond the training cohort.
- **Incremental training** — the `update` subcommand appends new labelled samples and retrains from the full feature cache, so deployed models stay current as cohorts grow.
- **Contributor guardrails** — [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md) documents quality standards, C++ style rules, and extension patterns for both human and AI contributors. The installation process also autogenerates type stubs, so the package comes with first-class type hinting support when used in Python.

---

## Quick start

Start by creating a manifest CSV that lists your samples:

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

See the [Usage](usage.md) page for full details.

Three feature classes drive the classifier:

| Feature class | Key signal |
|---|---|
| **End-motif frequencies** | Nuclease preference (DNase I / CAD balance shifts in disease) |
| **Fragment length** | Nucleosomal occupancy and chromatin remodelling |
| **CpG methylation** | Differential methylation at tissue-specific loci |

<!-- TODO this probably needs a citation - is the Li paper sufficient that we have in scratch/? -->
cfDNA fragments are generated by nuclease cleavage; different nucleases prefer different DNA sequences at the cut site (the "end motif"). Fragment length reflects chromatin packaging: nucleosomes protect small stretches of DNA, creating characteristic peaks at ~167 bp (mono-nucleosomal) and ~280–400 bp (di-nucleosomal).

See the [Results](results.md) page for plots, per-sample feature data, and classification metrics from the chr21 test cohort.

