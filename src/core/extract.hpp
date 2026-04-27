#pragma once

#include <string>
#include <unordered_map>

#include <htslib/sam.h>

namespace extractors {

// a better implementation would
// store motif size with the count
// map directly
using MotifCounts = std::unordered_map<std::string, size_t>;

void
accumlate_motifs
(MotifCounts& counts, const bam1_t* b, size_t motif_sz);


struct CpGSiteCount {
  size_t methylated = 0;
  size_t unmethylated = 0;
};
using CpGMap = std::unordered_map<int64_t, CpGSiteCount>;

void
accumulate_methylation
(CpGMap& cpg_sites, const bam1_t* b);

}

