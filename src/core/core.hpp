#pragma once

#include <cstddef>
#include <filesystem>
#include <string_view>
#include <unordered_map>

#include "stats.hpp"

namespace cfextract {

struct RegionMetrics {
  // end motifs: k-mer string -> read count
  std::unordered_map<std::string, size_t> end_motifs;
  FragLenBins frag_len_hist{};
  // summary statistics derived from frag_len_hist (NaN if no fragments observed)
  FragLenStats fl_stats;
  // CpG methylation summary statistics (NaN if no covered sites)
  MethylationStats methylation;
  // read start/end position histograms: (pos / 100_000) -> count, for summary plots only
  std::unordered_map<int64_t, size_t> start_pos_hist;
  std::unordered_map<int64_t, size_t> end_pos_hist;
};


RegionMetrics
extract_metrics
(
  std::filesystem::path aln_fp,
  std::string_view contig,
  int64_t start = 0,
  int64_t stop = -1,
  size_t motif_sz = 4
);

}
