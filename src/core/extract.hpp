#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string_view>
#include <unordered_map>

namespace cfextract {

struct CpGSiteCount {
  size_t methylated   = 0;
  size_t unmethylated = 0;
};

struct MethylationStats {
  double mean       = std::numeric_limits<double>::quiet_NaN();
  double entropy    = std::numeric_limits<double>::quiet_NaN();
  double frac_high  = std::numeric_limits<double>::quiet_NaN();
  double frac_low   = std::numeric_limits<double>::quiet_NaN();
  double median     = std::numeric_limits<double>::quiet_NaN();
};

struct FragLenStats {
  double mean                = std::numeric_limits<double>::quiet_NaN();
  double std_dev             = std::numeric_limits<double>::quiet_NaN();
  double frac_subnucleosomal = std::numeric_limits<double>::quiet_NaN();
  double frac_nucleosomal    = std::numeric_limits<double>::quiet_NaN();
  double ratio_mono_di       = std::numeric_limits<double>::quiet_NaN();
};

struct RegionMetrics {
  // end motifs: 4-mer string → read count
  std::unordered_map<std::string, size_t> end_motifs;
  // fragment lengths: index = length in bp, value = count; >1000 bp clamped to index 1000
  std::array<size_t, 1001> frag_len_hist{};
  // summary statistics derived from frag_len_hist (NaN if no fragments observed)
  FragLenStats fl_stats;
  // CpG methylation summary statistics (NaN if no covered sites)
  MethylationStats methylation;
  // read start/end position histograms: (pos / 100_000) → count, for summary plots only
  std::unordered_map<int64_t, size_t> start_pos_hist;
  std::unordered_map<int64_t, size_t> end_pos_hist;
};


RegionMetrics
extract_metrics
(
  std::filesystem::path aln_fp,
  std::string_view contig,
  int64_t start = 0,
  int64_t stop  = -1
);

}
