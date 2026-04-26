#include "extract.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "access.hpp"
#include "features.hpp"
#include "stats.hpp"


namespace {

void
accumulate_methylation
(const bam1_t* b, std::unordered_map<int64_t, cfextract::CpGSiteCount>& cpg_sites)
{
  const uint8_t* xm_ptr = bam_aux_get (b, "XM");
  if (!xm_ptr || (*xm_ptr != static_cast<uint8_t> ('Z'))) {
    return;
  }
  const char* xm = bam_aux2Z (xm_ptr);
  int64_t ref_p = b->core.pos;
  const uint32_t* cig = bam_get_cigar (b);
  size_t xi = 0;
  for (uint32_t ci = 0; ci < static_cast<uint32_t> (b->core.n_cigar); ++ci) {
    const uint32_t op = bam_cigar_op (cig[ci]);
    const uint32_t len = bam_cigar_oplen (cig[ci]);
    const auto type = bam_cigar_type (op);
    if ((type & 1) && (type & 2)) {
      for (uint32_t k = 0; k < len; ++k, ++xi, ++ref_p) {
        if (xm[xi] == 'Z') {
          ++cpg_sites[ref_p].methylated;
        }
        else if (xm[xi] == 'z') {
          ++cpg_sites[ref_p].unmethylated;
        }
      }
    }
    else if (type & 1) {
      xi += static_cast<size_t> (len);
    }
    else if (type & 2) {
      ref_p += static_cast<int64_t> (len);
    }
  }
}

}


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

  auto aln = ha::open_aln (aln_fp);

  const int32_t tid = sam_hdr_name2tid (aln.hdr, std::string (contig).c_str());
  if (tid < 0) {
    throw std::runtime_error ("unknown contig: " + std::string (contig));
  }
  if (stop < 0) {
    stop = static_cast<int64_t> (aln.hdr->target_len[static_cast<size_t> (tid)]);
  }

  std::unique_ptr<hts_itr_t, void(*)(hts_itr_t*)> reg_iter {
    ha::get_region (aln, {tid, start, stop}), hts_itr_destroy
  };
  std::unique_ptr<bam1_t, void(*)(bam1_t*)> b {bam_init1(), bam_destroy1};

  constexpr size_t motif_sz = 4;
  constexpr int64_t pos_bin_sz = 100000;

  RegionMetrics out;
  std::unordered_map<int64_t, CpGSiteCount> cpg_sites;

  while (true) {
    const auto itr_rc = sam_itr_next (aln.f, reg_iter.get(), b.get());
    // sam_itr_next returns -1 for EOF and < -1 for errors; both are terminal.
    if (itr_rc < 0) {
      break;
    }
    if ((b->core.flag & 3852) || (b->core.flag & 3) != 3) {
      continue;
    }

    // end motif
    ++out.end_motifs[fx::get_end_motif (b.get(), motif_sz)];

    // fragment length — leftmost read of pair only (isize > 0 avoids double-counting).
    // isize (TLEN) is set by the aligner; a more robust approach would derive fragment
    // length from PNEXT + cigar_ref_length(MC tag), which is independent of aligner
    // correctness. For Bismark PE output TLEN is reliable, but this is a known limitation.
    if (b->core.isize > 0) {
      const auto fl = static_cast<size_t> (b->core.isize);
      out.frag_len_hist[std::min (fl, size_t{1000})]++;
    }

    // methylation from XM Bismark tag — accumulated per-site; stats computed post-loop
    accumulate_methylation (b.get(), cpg_sites);

    // position distributions (100 kbp bins, summary visualization only)
    ++out.start_pos_hist[b->core.pos / pos_bin_sz];
    ++out.end_pos_hist[bam_endpos (b.get()) / pos_bin_sz];
  }

  // cpg_sites is local and freed here; only the 5 summary scalars cross the C++/Python boundary
  out.methylation = compute_meth_stats (cpg_sites);
  out.fl_stats    = compute_fl_stats (out.frag_len_hist);

  return out;
}

}
