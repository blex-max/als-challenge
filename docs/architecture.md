# Architecture

The pipeline has two layers connected by the `RegionMetrics` boundary type. Everything upstream of the boundary is C++; everything downstream is pure Python.

```
cfextract (C++ / pybind11)
│  src/core/core.{cpp,hpp}      — RegionMetrics struct; extract_metrics() entry point
│  src/core/access.{cpp,hpp}    — safe BAM/index file handle
│  src/core/stats.{cpp,hpp}     — binned data -> summary stats
│  src/core/extract.{cpp,hpp}   — per-read motif and methylation accumulators
│  src/bindings/python/         — bindings for python usage
│
│  ── RegionMetrics boundary ──────────────────────────────────────────────────
│
cfclassify (Python)
   cfclassify/types.py           — Data schemas (Features, Sample, ModelBundle)
   cfclassify/features.py        — metrics_to_features()
   cfclassify/classify.py        — LOO-CV, final model training, inference
   cfclassify/model.py           — model bundle and feature cache I/O
   cfclassify/plots.py           — matplotlib figures
   cfclassify/__main__.py        — CLI entry point (train/predict/update)
```

The boundary is deliberately narrow. `RegionMetrics` carries only compact, pre-reduced data: histograms (fixed-size arrays and small maps) and 5-element stats structs. Raw read data does not cross the C++/Python interface.

Python bindings are provided such that the extracted features can be used with common frameworks for downstream analysis, but the C++ core could equally be driven from R, Julia, or any other language with a C FFI. Additionally, when processing whole-genome BAM, I/O throughput will likely dominate total runtime. Processing alignment data directly in C++ with htslib is significantly more performant than using e.g. pysam due to inevitable Python overhead per record. This has not been explicitly benchmarked but the speedup is likely in the tens with respect to order of magnitude.

---

## Core Extractor Scalability

### Time complexity

Per `extract_metrics()` call: **O(R × L)** where R = reads in the region, L = read length. The inner per-read cost is O(1) for end-motif lookup and O(L) for CIGAR walking to parse the Bismark `XM` methylation tag. Post extraction computation is O(N log N) for the median computation over N covered CpG sites and O(1001) for fragment-length statistics. These latter demands are negligible relative to the BAM streaming loop.

### Memory model

**Peak runtime memory (inside `extract_metrics()`):**

The `cpg_sites` accumulation map is the dominant runtime allocation. The data stored for a given CpG site requires 24 bytes. Accounting for the overhead of the map itself, that figure grows to approximately 56 bytes per CpG. As such the expected memory usage for e.g. ~5 million covered CpG sites is 5M × 56 ~= 280 MB. On the subsampled chr12 Li data used for building this pipeline, memory usage peaks at about 140 KB for a typical chr21 alignment. This memory is freed before the function returns, the Python side is never required to handle it.

**Memory at the C++/Python boundary (returned `RegionMetrics`):**

| Structure | Size |
|---|---|
| `end_motifs` | ≤ 256 entries (all 4<sup>4</sup> 4-mers), ~8 KB |
| `frag_len_hist` | 1001 × 8 B = ~8 KB fixed |
| `fl_stats` | 5 doubles = 40 B |
| `methylation` | 5 doubles = 40 B |
| `start_pos_hist` / `end_pos_hist` | O(contig_len / 100 kbp): ~4 KB chr21 |

Python receives only the compact `RegionMetrics` struct regardless of coverage depth, i.e. only ~20 KB is serialised across the C++/Python boundary per call.

### Measured performance

End-to-end `cfextract.extract_features()` across 22 chr21 BAMs on an Apple M-series machine and using a Release build of this pipeline gave the following performance results:

| Metric | Mean ± std | Median |
|---|---|---|
| Wall time | 0.10 ± 0.00 s | 0.09 s |
| Peak RSS increase | 5.9 ± 0.3 MB | 5.8 MB |

Reproduce with `python scripts/bench_extract.py --repeat 3`.

---

## Extending the tool

Adding a new feature to the pipeline requires changes in exactly two places:

1. **C++ side** (`src/core/core.hpp` and `src/core/core.cpp`): add the new field to `RegionMetrics` and compute its final value inside `extract_metrics()`. If the feature requires per-read accumulation, add a helper in `src/core/extract.{hpp,cpp}` (alongside `accumlate_motifs` and `accumulate_methylation`) and call it from the per-read loop in `extract_metrics()`.

2. **Python side** (`cfclassify/features.py`, `metrics_to_features()`): read the new attribute from the `RegionMetrics` object and write derived value(s) into the `Features` dict. If the result is a flat numeric scalar (like `fl_mean` or `methylation_entropy`), `build_feature_matrix()` in `classify.py` picks it up automatically — no further changes needed. If it is a histogram or dict (like `end_motifs`), add explicit handling to `build_feature_matrix()` as well.

No changes are needed to the classification, model I/O, or plotting infrastructure.

---

## Engineering limitations

<!-- TODO separate the classifer limiations (bullet point 3) out - that should be in the classifier section it is not an eng concern -->
- **Model storage**: `ModelBundle` is serialised with `joblib` (pickle). There is no versioning, schema migration, or compatibility guarantee across Python or sklearn releases. A bundle saved with one sklearn version may not load cleanly after a major upgrade.
- **Model scalability**: the extractor scales to large cohorts (I/O-bound, memory bounded by region size). The classifier does not: `LogisticRegression` has no `partial_fit`, so `update` always retrains from scratch on all cached feature vectors. This is fast at cfDNA cohort sizes (seconds on CPU) but would not scale to thousands of samples without switching to a solver that supports incremental updates.
- **Single-threaded extraction**: `extract_metrics()` processes one BAM serially. The CLI uses `joblib` to parallelise across samples at the process level, but per-BAM throughput is single-threaded. This is sufficient for chr21 BAMs (~0.1 s each) but would be a bottleneck for whole-genome runs without introducing intra-file parallelism.
- **Single target region**: The current implementation of feature extraction takes a single target region, and the Python CLI tool simply defaults to targeting chr21. A python side programmer would need to make multiple calls `extract_metrics` to target multiple regions or the entire genome. A future implementation could efficiently accumulate data from multiple target regions using existing htslib multi-region iterator infrastructure.
