#include "features.hpp"


namespace feature_extractors {

std::string
get_end_motif
(const bam1_t* b, size_t motif_sz)
{
  std::string motif;

  const auto bseq = bam_get_seq (b);

  if (b->core.flag & BAM_FREVERSE) {
    // does the motif need to be reverse complemented?
    // is the reverse complement of the motif to be treated as degenerate?
    const auto seq_len = b->core.l_qseq;
    for (size_t i = seq_len - motif_sz; i < seq_len; ++i) {
      motif.push_back (seq_nt16_str[bam_seqi(bseq, i)]);
    }
  }
  else {
    for (size_t i = 0; i < motif_sz; ++i) {
      motif.push_back (seq_nt16_str[bam_seqi(bseq, i)]);
    }
  }

  return motif;
}

}
