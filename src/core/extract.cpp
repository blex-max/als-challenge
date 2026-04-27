#include "extract.hpp"

#include <algorithm>

#include "util.hpp"


namespace extractors {

void
accumlate_motifs
(MotifCounts& counts, const bam1_t* b, size_t motif_sz)
{
  std::string motif;

  const auto bseq = bam_get_seq (b);
  const auto seq_len = static_cast<size_t> (b->core.l_qseq);

  if (b->core.flag & BAM_FREVERSE) {
    for (size_t i = seq_len - motif_sz; i < seq_len; ++i) {
      motif.push_back (seq_nt16_str[bam_seqi(bseq, i)]);
    }
    std::reverse (motif.begin(), motif.end());
    for (char& c : motif) {
      c = util::seq_complement (c);
    }
  }
  else {
    for (size_t i = 0; i < motif_sz; ++i) {
      motif.push_back (seq_nt16_str[bam_seqi(bseq, i)]);
    }
  }

  ++counts[motif];
}

void
accumulate_methylation
(CpGMap& cpg_sites, const bam1_t* b)
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

