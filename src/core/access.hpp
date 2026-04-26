#pragma once

#include <filesystem>

#include <htslib/hts.h>
#include <htslib/sam.h>

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

  AlnFile() = default;
  AlnFile(const AlnFile&) = delete;
  AlnFile& operator=(const AlnFile&) = delete;

  AlnFile (AlnFile&& o) noexcept {
    f = o.f;
    hdr = o.hdr;
    idx = o.idx;
    o.f = nullptr;
    o.hdr = nullptr;
    o.idx = nullptr;
  }

  ~AlnFile () {
    if (idx) {
      hts_idx_destroy (idx);
    }
    if (hdr) {
      sam_hdr_destroy (hdr);
    }
    if (f) {
      hts_close (f);
    }
  }
};

AlnFile
open_aln
(std::filesystem::path fp);

hts_itr_t*
get_region
(const AlnFile& aln, GenomicRegion reg);

};
