# ALS cfDNA Analysis

A cfDNA feature-extraction and exploratory classification pipeline for distinguishing ALS and control samples in a small research cohort using bisulfite-sequenced plasma cell-free DNA. Cell-free DNA (cfDNA) in blood plasma is shed by apoptotic and necrotic cells throughout the body. Because cfDNA carries epigenetic and fragmentation signatures from its tissue of origin, it is being investigated as a liquid-biopsy substrate for neurodegenerative disease, including ALS.

---

## Pipeline at a glance

A high-performance C++ core (`cfextract`) handles BAM I/O via htslib and accumulates per-read statistics (end motifs, fragment lengths, and CpG methylation) into a compact `RegionMetrics` boundary object. A Python modelling layer using sklearn (`cfclassify`) converts that boundary object to a flat feature vector and trains or applies an L2 logistic regression classifier. On the chr21 cohort of 22 bisulfite-sequenced plasma BAMs (12 ALS, 10 CTRL), the pipeline extracts 266 features per sample and reaches 0.64 accuracy under LOO-CV.

### Engineering Highlights

- **CI on every push** — a series of code analysis tools, and unit tests, run on every commit; failures surface in GitHub Checks, so regressions and bugs are caught before they reach a production environment. Analysis includes code complexity/maintainability assessment.
- **Portable container** — pre-built image available via package registry (`docker pull ghcr.io/blex-max/als-challenge:latest`) minimises end-user issues; the full pipeline should run reproducibly on any machine with a single pull.
- **Language-agnostic extraction core** — core feature extraction and file I/O are handled in performant cpp. Python bindings are provided out of the box, but analysts can drive analysis from Python, R, Julia, or any language for which bindings can be made, without having to reimplement key functionality.
- **Memory scales with region size only, not read depth** — peak footprint is bounded by the number of unique CpG sites in the target region, not read count; a deeply sequenced chr21 BAM uses the same memory as a shallow one, keeping extraction tractable on standard hardware.
- **Fast extraction** — `cfextract.extract_features()` completes chr21 extraction in **0.10 ± 0.00 s** wall time, 5.9 ± 0.3 MB peak memory usage across 22 chr21 BAMs.
- **Multiple Entrypoints** — the `train` subcommand fits and persists a final model bundle; `predict` applies it to new samples without retraining, enabling use and testing beyond the training cohort.
- **Incremental training** — the `update` subcommand appends new labelled samples and retrains from the full feature cache, so deployed models stay current as cohorts grow.
- **Contributor guardrails** — [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md) documents quality standards, C++ style rules, and extension patterns for both human and AI contributors. The installation process also autogenerates type stubs, so the package comes with first-class type hinting support when used in Python.

See the [Architechture](architecture.md) page for a more in-depth discussion of the design.

---

## Quick start

To train the model, start by creating a manifest CSV that lists your samples:

```csv
sample_id,bam_path,label
ALS_001,bams/ALS_001.bam,als
CTRL_001,bams/CTRL_001.bam,ctrl
```

Simple usage is then as follows:

**Option A — Docker (no local build required)**

```bash
docker pull ghcr.io/blex-max/als-challenge:latest

# Train: runs LOO-CV evaluation and saves a deployable model to /results/model.pkl
docker run --rm \
  -v /path/to/bams:/data/bams \
  -v /path/to/samples.csv:/data/samples.csv \
  -v /path/to/results:/results \
  ghcr.io/blex-max/als-challenge:latest \
  train --manifest path/to/manifest.csv --out-dir /results

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

cfclassify train --manifest path/to/manifest.csv --out-dir results/
cfclassify predict --bam new_patient.bam --model-path results/model.pkl --sample-id new_patient
```

See the [Usage](usage.md) page for full details.

Three feature classes drive the classifier:

| Feature class | Key signal |
|---|---|
| **End-motif frequencies** | Differential appearance of certain k-mers at fragment ends |
| **Fragment length** | Length distribution of cfDNA fragments |
| **CpG methylation** | Differential methylation between samples |

See the [Results](results.md) page for plots, per-sample feature data, and classification metrics from the chr21 test cohort.


---

## References

1. Caggiano C, Celona B, Garton F, et al. *Comprehensive cell type decomposition of circulating cell-free DNA with CelFiE*. Nature Communications. 2021;12:2717. https://doi.org/10.1038/s41467-021-22901-x

2. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. *Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin*. Cell. 2016;164(1-2):57–68. https://doi.org/10.1016/j.cell.2015.11.050

3. Ding SC, Lo YMD. *Cell-Free DNA Fragmentomics in Liquid Biopsy*. Diagnostics. 2022;12(4):978. https://doi.org/10.3390/diagnostics12040978

4. Moss J, Magenheim J, Neiman D, et al. *Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease*. Nature Communications. 2018;9:5068. https://doi.org/10.1038/s41467-018-07466-6


