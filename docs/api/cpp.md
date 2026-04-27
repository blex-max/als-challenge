# C++ API — cfextract

The C++ core lives in `src/core/`. It is exposed to Python via pybind11 bindings in `src/bindings/python/bindings.cpp`.

!!! note "Coordinate convention"
    All genomic coordinates are **0-based, half-open** (`[start, stop)`), matching the SAM/BAM convention used by htslib. `stop = -1` means "to the end of the contig".

---

## Primary entry point

```cpp
namespace cfextract {

struct MethylationStats {
  double mean, entropy, frac_high, frac_low, median;  // all NaN if no covered sites
};

struct FragLenStats {
  double mean, std_dev, frac_subnucleosomal, frac_nucleosomal, ratio_mono_di;  // all NaN if no fragments
};

struct RegionMetrics {
  std::unordered_map<std::string, size_t> end_motifs;   // k-mer → read count
  std::array<size_t, 1001>               frag_len_hist; // index = length in bp; ≥1000 clamped to 1000
  FragLenStats                           fl_stats;
  MethylationStats                       methylation;
  std::unordered_map<int64_t, size_t>   start_pos_hist; // (pos / 100_000) → count
  std::unordered_map<int64_t, size_t>   end_pos_hist;
};

RegionMetrics extract_metrics(
  std::filesystem::path aln_fp,
  std::string_view      contig,
  int64_t               start    = 0,
  int64_t               stop     = -1,
  size_t                motif_sz = 4
);

} // namespace cfextract
```

### `cfextract::extract_metrics()`

Opens the BAM at `aln_fp`, iterates all properly-paired primary reads in `[start, stop)` on `contig`, and returns a fully-populated `RegionMetrics`. The function:

1. Skips unmapped, secondary, QC-failed, and duplicate reads.
2. Only processes reads where both mates are mapped (`FLAG & 0x3 == 0x3`).
3. Extracts the k-mer end motif from each read (reverse-complemented for reverse-strand reads).
4. Records fragment length from `TLEN` for the leftmost read of each pair (`TLEN > 0`); lengths ≥ 1000 bp are clamped to index 1000 in `frag_len_hist`.
5. Accumulates per-position CpG counts from the Bismark `XM` auxiliary tag.
6. Bins read start and end positions into 100 kbp bins for position histograms.
7. Computes `FragLenStats` and `MethylationStats` from the accumulated data before freeing the internal `cpg_sites` map.

Python signature (via pybind11):

```python
cfextract.extract_features(
    aln_filepath: str,
    contig: str,
    region_start: int = 0,
    region_end: int = -1,
    motif_sz: int = 4,
) -> cfextract.RegionMetrics
```

---

## BAM access — `htsacc`

```cpp
namespace htsacc {

struct AlnFile {                        // non-copyable RAII wrapper
  htsFile*   f   = nullptr;
  sam_hdr_t* hdr = nullptr;
  hts_idx_t* idx = nullptr;
  ~AlnFile();                           // frees all three handles in reverse-acquisition order
};

AlnFile   open_aln  (std::filesystem::path fp);
hts_itr_t* get_region(const AlnFile& aln, GenomicRegion reg);

} // namespace htsacc
```

`AlnFile` wraps the three htslib handles required for indexed BAM access. The destructor frees all three in reverse-acquisition order.

`open_aln()` opens the file and loads the BAI/CSI index, throwing on any failure. `get_region()` returns an `hts_itr_t*` iterator for the specified region; the caller is responsible for freeing it with `hts_itr_destroy()`.

---

## End-motif extraction — `extractors`

```cpp
namespace extractors {

using MotifCounts = std::unordered_map<std::string, size_t>;

void accumlate_motifs(MotifCounts& counts, const bam1_t* b, size_t motif_sz);

} // namespace extractors
```

`accumlate_motifs()` updates a running k-mer tally with the end motif from a single read:

- **Forward-strand reads** (`BAM_FREVERSE` not set): first `motif_sz` bases of the read sequence.
- **Reverse-strand reads** (`BAM_FREVERSE` set): last `motif_sz` bases, reversed and complemented — this gives the 5′ end of the complementary strand, i.e. the actual fragment end.

The accumulated `MotifCounts` map is transferred directly into `RegionMetrics::end_motifs` by `extract_metrics()`.

Bismark bisulfite alignment converts unmethylated C to T, so motif distributions are T/A-biased compared to non-bisulfite cfDNA.

---

## Summary statistics — `cfextract`

```cpp
namespace cfextract {

MethylationStats compute_meth_stats(const extractors::CpGMap& sites);
FragLenStats     compute_fl_stats  (const FragLenBins& hist);

} // namespace cfextract
```

Both functions return all-NaN structs if their inputs are empty (no fragments or no covered CpG sites). The NaN sentinel propagates to Python via the pybind11 bindings and is detected in `metrics_to_features()` — absent data becomes an absent key rather than a NaN value in the `Features` dict.

### `compute_fl_stats(hist)`

Computed from the 1001-bin length histogram in a single pass:

- **mean**, **std_dev**: standard online formulas — `Σ(i × hist[i]) / N` and `√(Σ(hist[i] × (i − mean)²) / N)`.
- **frac_subnucleosomal**: `sum(hist[0..119]) / N` — fragments shorter than 120 bp.
- **frac_nucleosomal**: `sum(hist[120..200]) / N` — mono-nucleosomal peak.
- **ratio_mono_di**: `frac_nucleosomal / frac_dinucleosomal`, where di-nucleosomal = `sum(hist[280..400]) / N`. This ratio is the strongest single discriminator in the chr21 cohort (effect size 0.84).

### `compute_meth_stats(sites)`

Computed from the `extractors::CpGMap` (`unordered_map<int64_t, CpGSiteCount>`) internal accumulator:

- Per-site methylation rate: `methylated / (methylated + unmethylated)`; sites with zero coverage are skipped.
- **mean**: coverage-weighted — `Σ(methylated) / Σ(total)` — equivalent to treating all reads equally rather than all sites equally.
- **entropy**: mean binary entropy `H(r) = −r log₂r − (1−r) log₂(1−r)` across sites; boundary values (r = 0 or 1) treated as 0.
- **frac_high** / **frac_low**: fraction of sites with rate > 0.8 / < 0.1.
- **median**: 50th percentile of sorted per-site rates.
