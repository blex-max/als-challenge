#include "features.hpp"

#include <algorithm>


namespace {

char
complement (char c) noexcept
{
  switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default:  return 'N';
  }
}

}


namespace feature_extractors {

std::string
get_end_motif
(const bam1_t* b, size_t motif_sz)
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
      c = complement (c);
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
