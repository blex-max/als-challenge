#include "access.hpp"

#include <stdexcept>

#include <htslib/hts.h>
#include <htslib/sam.h>

namespace htsacc {

AlnFile
open_aln
(std::filesystem::path fp)
{
  AlnFile aln;

  aln.f = hts_open (fp.c_str(), "r");
  if (!aln.f) {
    throw std::runtime_error ("failed to open: " + fp.string());
  }

  aln.hdr = sam_hdr_read (aln.f);
  if (!aln.hdr) {
    throw std::runtime_error ("failed to read header: " + fp.string());
  }

  aln.idx = sam_index_load (aln.f, fp.c_str());
  if (!aln.idx) {
    throw std::runtime_error ("failed to load index: " + fp.string() + " (ensure .bai/.csi exists)");
  }

  return aln;
}

hts_itr_t*
get_region
(const AlnFile& aln, GenomicRegion reg)
{
  return sam_itr_queryi (aln.idx, reg.tid, reg.start, reg.stop);
}

}
