# ALS cfDNA Analysis

A cfDNA feature-extraction and classification pipeline for distinguishing ALS patients from healthy controls using bisulfite-sequenced plasma cell-free DNA.

---

## Background

Cell-free DNA (cfDNA) in blood plasma is shed by apoptotic and necrotic cells throughout the body. Because cfDNA carries epigenetic and fragmentation signatures from its tissue of origin, it can serve as a liquid biopsy for a range of conditions — including neurodegenerative diseases such as ALS.

The question is: can we distinguish ALS patients from healthy controls using features extracted from bisulfite-sequenced cfDNA BAM files?

The reference cohort comes from Li et al. (2021) (PRJNA691320): 22 samples (12 ALS, 10 CTRL), Illumina NovaSeq, Bismark-aligned, chr21 only for this analysis.

---

## Pipeline at a glance

```
BAM files (12 ALS + 10 CTRL, chr21, bisulfite-sequenced)
         │
         ▼
  cfextract  (C++ / pybind11 / htslib)
  ─────────────────────────────────────
  • Per-fragment: end 4-mer motif, length, XM methylation tag
  • Histograms accumulated in O(bins) memory — not O(reads)
  • Returns compact RegionMetrics struct
         │
         ▼
  cfanalysis  (pure Python)
  ─────────────────────────
  • Normalises end-motif counts → frequencies
  • Extracts scalar summary statistics (fragment length, methylation)
  • Runs LOO-CV with L2 logistic regression + inner GridSearchCV
  • Produces plots and classification report
```

---

## Engineering highlights

- **CI on every push** — A series of code analysis tools, and unit tests, run on every commit; failures surface in GitHub Checks, so regressions and bugs are caught before they reach a production environment.
- **Portable container** — pre-built image available via package registry (`docker pull ghcr.io/blex-max/als-challenge:latest`) eliminates end-user build issues; the full pipeline runs reproducibly on any machine with a single pull.
- **Language-agnostic extraction core** — Core feature extraction and file I/O are handled in performant cpp. Python bindings are provided out of the box, but analysts can drive analysis from Python, R, Julia, or any language for which bindings can be made, without having to reimplement key functionality.
- **Memory scales with region size only, not read depth** — peak footprint is bounded by the number of unique CpG sites in the target region, not read count; a deeply sequenced chr21 BAM uses the same memory as a shallow one, keeping extraction tractable on standard hardware.
- **Fast extraction with benchmarks** — fragment stats complete in ~2 µs, methylation at 50k sites in ~513 µs; benchmarks are embedded in the test suite and are verifiable on any machine.
- **Contributor guardrails** — [CONTRIBUTING.md](https://github.com/blex-max/als-challenge/blob/main/CONTRIBUTING.md) documents quality standards, C++ style rules, and extension patterns for human and AI contributors alike

---

Three feature classes drive the classifier:

| Feature class | Key signal | Chr21 performance |
|---|---|---|
| **End-motif frequencies** | Nuclease preference (DNase I / CAD balance shifts in disease) | Moderate |
| **Fragment length** | Nucleosomal occupancy and chromatin remodelling | **Strong** (`fl_ratio_mono_di` effect size 0.89) |
| **CpG methylation** | Differential methylation at tissue-specific loci | Weak (chr21 is constitutively methylated) |

---

## Current results (chr21 cohort)

LOO-CV with L2 logistic regression on 22 samples:

| Metric | All features | Fragment-only |
|---|---|---|
| Accuracy | 0.64 | 0.68 |
| Macro F1 | 0.62 | 0.66 |

Fragment length features — especially `fl_ratio_mono_di` (mono-nucleosomal / di-nucleosomal count ratio) — dominate the signal. Methylation features add noise on chr21 because this chromosome is dominated by satellite repeats that are essentially fully methylated in all samples. See [Classification](classification.md) for the full feature analysis, classifier justification, and LOO-CV methodology.

!!! note "Full-genome projection"
    Chr21 is ~1.5% of the autosomal genome — extending the analysis to the full genome is expected to substantially increase classification power. Notably, methylation features are uninformative here but likley carry genuine signal across the rest of the genome.

---

## Quick start

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .

cfanalysis \
  --manifest data/samples.csv \
  --output-dir results/
```

Or via Docker:

```bash
docker pull ghcr.io/blex-max/als-challenge:latest
docker run --rm \
  -v /path/to/bams:/data/bams \
  -v /path/to/samples.csv:/data/samples.csv \
  -v /path/to/results:/results \
  ghcr.io/blex-max/als-challenge:latest \
  --manifest /data/samples.csv --output-dir /results
```

See the [Usage](usage.md) page for full installation and options.
