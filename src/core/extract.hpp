#pragma once

#include <filesystem>
#include <string_view>
#include <unordered_map>

namespace cfextract {

// columnar table representation
// provides best memory layout
// for most commonly expected access
// pattern - calculation per feature
// (rather than per fragment).
// struct FragmentMetrics {
//   std::vector<std::string> qname;
//   std::vector<int64_t> frag_start;
//   std::vector<int64_t> frag_end;
//   std::vector<int32_t> frag_len;
//   std::vector<std::array<char, 4>> start_4mer;
//   std::vector<std::array<char, 4>> end_4mer;
// };

// struct CpGMetrics {
//   std::vector<meth>
// }

// struct SampleMetrics {
//   std::string sample_id;
//   std::string bam_path;
//   /* more sample level info */

//   std::vector<FragmentMetrics> region_metrics;
//   std::vector<std::string> regions;  // contig info, perhaps a struct
// };


// struct CpGSiteCount {
//   size_t methylated;
//   size_t unmethylated;
// };

struct RegionMetrics {
  std::unordered_map<std::string, size_t> end_motifs;  // count of each motif seen
  // std::unordered_map<int64_t, CpGSiteCount> cpg_counts;  // keyed by genomic position
};


RegionMetrics
extract_metrics
(
  std::filesystem::path aln_fp,
  std::string_view contig,
  int64_t start = 0,
  int64_t stop  = -1   // -1 -> whole contig
  /* opts */
);

}

