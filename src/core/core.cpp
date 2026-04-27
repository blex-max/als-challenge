#include "core.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "access.hpp"
#include "stats.hpp"

constexpr int64_t POS_BIN_SZ = 100000;  // hardcoded for simplicity

namespace cfextract {

namespace ha = htsacc;
namespace ext = extractors;

RegionMetrics
extract_metrics
(
  std::filesystem::path aln_fp,
  std::string_view contig,
  int64_t start,
  int64_t stop,
  size_t motif_sz
)
{
  auto aln = ha::open_aln (aln_fp);

  const int32_t tid = sam_hdr_name2tid (aln.hdr, std::string (contig).c_str());
  if (tid < 0) {
    throw std::runtime_error ("unknown contig: " + std::string (contig));
  }
  if (stop < 0) {
    stop = static_cast<int64_t> (aln.hdr->target_len[static_cast<size_t> (tid)]);
  }

  const auto reg_iter = ha::get_region (aln, {tid, start, stop});
  auto b = bam_init1();


  RegionMetrics out;
  ext::CpGMap cpg_sites;

  while (true) {
    const auto itr_rc = sam_itr_next (aln.f, reg_iter, b);
    // sam_itr_next returns -1 for EOF and < -1 for errors.
    // In an extended implementation proper
    // error state handling would be nice.
    if (itr_rc < 0) {
      break;
    }
    if ((b->core.flag & 3852) || (b->core.flag & 3) != 3) {
      continue;
    }

    // end motif
    ext::accumlate_motifs (out.end_motifs, b, motif_sz);

    // fragment length — leftmost read of pair only (isize > 0 avoids double-counting).
    // isize (TLEN) is set by the aligner; a more robust approach would derive fragment
    // length from PNEXT + cigar_ref_length(MC tag). For Bismark PE output TLEN is
    // ostensibly reliable, but this is a known limitation.
    if (b->core.isize > 0) {
      const auto fl = static_cast<size_t> (b->core.isize);
      out.frag_len_hist[std::min (fl, size_t{1000})]++;
    }

    // methylation from XM Bismark tag — accumulated per-site; stats computed post-loop
    ext::accumulate_methylation (cpg_sites, b);

    ++out.start_pos_hist[b->core.pos / POS_BIN_SZ];
    ++out.end_pos_hist[bam_endpos (b) / POS_BIN_SZ];
  }

  out.methylation = compute_meth_stats (cpg_sites);
  out.fl_stats = compute_fl_stats (out.frag_len_hist);

  hts_itr_destroy (reg_iter);
  bam_destroy1 (b);

  return out;
}

}
