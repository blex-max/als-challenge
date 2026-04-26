# Architecture

The pipeline has two layers connected by the `RegionMetrics` boundary type. Everything upstream of the boundary is C++; everything downstream is pure Python.

```
cfextract   (C++ / pybind11)
│  src/core/access.{cpp,hpp}    — safe BAM/index file handle
│  src/core/features.{cpp,hpp}  — per-read end-motif extraction
│  src/core/stats.{cpp,hpp}     — histogram → scalar summary stats
│  src/core/extract.{cpp,hpp}   — top-level extraction loop → RegionMetrics
│  src/bindings/python/         — bindings for python usage
│
│  ── RegionMetrics boundary ──────────────────────────────────────────────────
│
cfanalysis   (pure Python)
   cfanalysis/types.py           — Data schemas
   cfanalysis/features.py        — metrics_to_features()
   cfanalysis/classify.py        — LOO-CV with inner GridSearchCV
   cfanalysis/plots.py           — matplotlib figures
   cfanalysis/__main__.py        — CLI entry point
```

The boundary is deliberately narrow. `RegionMetrics` carries only compact, pre-reduced data: histograms (fixed-size arrays and small maps) and 5-element stats structs. Raw read data does not cross the C++/Python interface.

Python bindings are provided such that the extracted features can be used with common frameworks for downstream analysis, but the C++ core can equally be driven from R, Julia, or any other language with a C FFI.

---

## Scalability

### Time complexity

Per `extract_metrics()` call: **O(R × L)** where R = reads in the region, L = mean read length. The inner per-read cost is O(1) for end-motif lookup and O(L) for CIGAR walking to parse the Bismark `XM` methylation tag. Post extraction computation is O(N log N) for the median computation over N covered CpG sites and O(1001) for fragment-length statistics. These latter demands are negligible relative to the BAM streaming loop.

### Space complexity

| Structure | Lifetime | Size |
|---|---|---|
| `end_motifs` | returned | ≤ 256 entries (all 4<sup>4</sup> 4-mers), ~8 KB |
| `frag_len_hist` | returned | 1001 × 8 B = ~8 KB fixed |
| `fl_stats` | returned | 5 doubles = 40 B |
| `methylation` | returned | 5 doubles = 40 B |
| `start_pos_hist` / `end_pos_hist` | returned | O(contig_len / 100 kbp): ~4 KB chr21, ~240 KB whole genome |
| `cpg_sites` (internal) | freed before return | O(covered CpG sites): ~140 KB chr21, up to ~280 MB whole genome |

`cpg_sites` is a local variable inside `extract_metrics()` and is freed before the function returns. Python receives only the 5-scalar `MethylationStats` struct (40 bytes) regardless of coverage depth or region size — the difference between serialising ~280 MB vs 40 bytes across the C++/Python boundary per call.

### Why C++ over pure Python

A pure-Python implementation using pysam with Welford-style streaming accumulators could compute the coverage-weighted methylation mean in O(1) space. However, four of the five methylation features (entropy, high/low fractions, median) require the **final per-site methylation rate** before they can be computed. A CpG position can be covered by many reads arriving non-consecutively in the BAM, so a site's rate is not known until all reads covering it have been processed.

The per-site accumulation map (`cpg_sites`) is therefore the minimal-state online algorithm for these features — it processes one read at a time and updates per-site counters, and cannot be replaced by running accumulators without approximating the median (e.g., Greenwald-Khanna sketch) and losing exact entropy. The map stays entirely in C++ and is discarded before returning; Python never touches it.

The secondary argument is raw BAM I/O throughput: htslib processes BAM records ~10–50× faster than pysam due to Python's GIL and ctypes overhead per record. For 50 GB BAMs at ~500 M reads, I/O dominates total runtime.

### Measured performance

Timings from the Catch2 benchmark suite (Apple M-series, Release build):

| Routine | Input | Time |
|---|---|---|
| `compute_fl_stats` | 1001-bin histogram | ~2 µs |
| `compute_meth_stats` | 50k CpG sites | ~513 µs |

Reproduce with `test-cfextract --benchmark-samples 100`.

---

## Extension point

Adding a new feature to the pipeline requires changes in only two places:

1. **Add a field to `RegionMetrics`** (C++, `src/core/extract.hpp`) and populate it in `src/core/extract.cpp`.
2. **Add corresponding entries in `cfanalysis/features.py::metrics_to_features()`** — flat numeric keys are automatically picked up by `build_feature_matrix()` in `classify.py`.

No changes are needed to the classification or plotting infrastructure.

```python title="cfanalysis/features.py — the extension point" linenums="1"
--8<-- "cfanalysis/features.py"
```

---

## Data flow

```
cfextract.extract_features(bam_path, contig)
    │  htslib BAM I/O
    │  per-read: end motif, fragment length, XM methylation accumulation
    │  on return: compute_fl_stats() + compute_meth_stats() + free cpg_sites
    ▼
RegionMetrics
    │
    ▼
metrics_to_features(metrics)          # cfanalysis/features.py
    │  normalise end-motif counts
    │  copy histogram + scalar stats
    │  NaN-check: absent data → absent key (not NaN value)
    ▼
Features dict  →  Sample(id, label, features)
    │
    ├─→  run_loo_cv(samples)           # cfanalysis/classify.py
    │        LOO outer loop
    │          ├─ _select_top_motifs() on training fold
    │          ├─ build_feature_matrix()
    │          ├─ StandardScaler.fit_transform() on training fold
    │          └─ GridSearchCV → predict held-out sample
    │        classification_report
    │
    └─→  plot_*(samples)               # cfanalysis/plots.py
             end motifs / frag lengths / methylation / position histograms
```
