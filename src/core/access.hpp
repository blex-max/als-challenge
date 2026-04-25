#pragma once

#include <filesystem>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string_view>

namespace htsacc {

struct GenomicRegion {
  int32_t tid;
  int64_t start;
  int64_t stop;
};

struct AlnFile {
  htsFile* f=nullptr;
  sam_hdr_t* hdr=nullptr;
  hts_idx_t* idx=nullptr;
};  // TODO destroy on close

AlnFile
open_aln
(std::filesystem::path fp);

hts_itr_t*
get_region
(const AlnFile& aln, GenomicRegion reg);

std::string_view
get_seq
(const bam1_t* b);

};
