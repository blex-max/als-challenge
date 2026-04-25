#include "extract.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "access.hpp"
#include "features.hpp"


namespace cfextract {

namespace ha = htsacc;

RegionMetrics
extract_metrics
(
  std::filesystem::path aln_fp,
  std::string_view contig,
  int64_t start,
  int64_t stop
)
{
  namespace fx = feature_extractors;

  RegionMetrics out;
  auto& motif_map = out.end_motifs;

  auto aln = ha::open_aln (aln_fp);

  const int32_t tid = sam_hdr_name2tid (aln.hdr, std::string (contig).c_str());
  if (tid < 0) {
    hts_idx_destroy (aln.idx);
    sam_hdr_destroy (aln.hdr);
    hts_close (aln.f);
    throw std::runtime_error ("unknown contig: " + std::string (contig));
  }
  if (stop < 0)
    stop = static_cast<int64_t> (aln.hdr->target_len[static_cast<size_t> (tid)]);

  auto reg_iter = ha::get_region (aln, {tid, start, stop});

  auto b = bam_init1();

  size_t motif_sz = 4;   // TODO propagate up as option to this function;

  // TODO multi region iter
  while (true) {
    const auto itr_rc = sam_itr_next (aln.f, reg_iter, b);
    if (itr_rc < -1) {
      // TODO error state
      break;
    }
    if (itr_rc < 0) {
      // region exhausted
      break;
    }
    if ((b->core.flag & 3852) || (b->core.flag & 3) != 3) {
      // bad read
      continue;
    }

    // extract features
    const auto end_motif = fx::get_end_motif(b, motif_sz);
    auto motif_bin = motif_map.find (end_motif);
    if (motif_bin == end (motif_map)) {
      motif_map[end_motif] = 1;
    }
    else {
      ++motif_map[end_motif];
    }
  }

  hts_itr_destroy (reg_iter);
  hts_idx_destroy (aln.idx);
  sam_hdr_destroy (aln.hdr);
  hts_close (aln.f);

  return out;
}

}
