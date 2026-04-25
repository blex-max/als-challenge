#include <catch2/catch_test_macros.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

#include <htslib/sam.h>

#include "core/access.hpp"
#include "core/extract.hpp"
#include "core/features.hpp"

// ---------------------------------------------------------------------------
// Helper: build a minimal synthetic bam1_t with a known sequence.
// data layout: 1-byte qname (null) + packed 4-bit sequence.
// bam_destroy1 calls free(b->data), so calloc ownership is correct.
// ---------------------------------------------------------------------------
static bam1_t*
make_read(const std::string& seq, bool reverse_strand)
{
  bam1_t* b       = bam_init1();
  b->core.l_qname = 1;
  b->core.n_cigar = 0;
  b->core.l_qseq  = static_cast<int32_t>(seq.size());
  b->core.flag    = reverse_strand
                      ? static_cast<uint16_t>(BAM_FREVERSE)
                      : uint16_t{0};

  const std::size_t data_sz = 1 + static_cast<std::size_t>((seq.size() + 1) / 2);
  b->data   = static_cast<uint8_t*>(std::calloc(data_sz, 1));
  b->l_data = static_cast<int>(data_sz);
  b->m_data = static_cast<uint32_t>(data_sz);

  uint8_t* packed = b->data + 1;
  for (std::size_t i = 0; i < seq.size(); ++i) {
    const uint8_t code = seq_nt16_table[static_cast<unsigned char>(seq[i])];
    if (i % 2 == 0) packed[i / 2]  = static_cast<uint8_t>(code << 4U);
    else             packed[i / 2] |= code;
  }
  return b;
}

// ---------------------------------------------------------------------------
// RegionMetrics
// ---------------------------------------------------------------------------

TEST_CASE("RegionMetrics default-constructed end_motifs is empty") {
  cfextract::RegionMetrics m;
  REQUIRE(m.end_motifs.empty());
}

// ---------------------------------------------------------------------------
// open_aln error path
// ---------------------------------------------------------------------------

TEST_CASE("open_aln throws std::runtime_error for non-existent file") {
  std::cerr << "[note] htslib error on the next line is expected\n";
  REQUIRE_THROWS_AS(
    htsacc::open_aln("/nonexistent/path/does-not-exist.bam"),
    std::runtime_error
  );
}

// ---------------------------------------------------------------------------
// get_end_motif
// ---------------------------------------------------------------------------

TEST_CASE("get_end_motif forward strand returns first N bases") {
  bam1_t* b = make_read("ACGTTTTT", false);
  REQUIRE(feature_extractors::get_end_motif(b, 4) == "ACGT");
  bam_destroy1(b);
}

TEST_CASE("get_end_motif reverse strand returns last N bases") {
  bam1_t* b = make_read("TTTTACGT", true);
  REQUIRE(feature_extractors::get_end_motif(b, 4) == "ACGT");
  bam_destroy1(b);
}

TEST_CASE("get_end_motif respects smaller motif size") {
  bam1_t* b = make_read("ACGTTTTT", false);
  REQUIRE(feature_extractors::get_end_motif(b, 2) == "AC");
  bam_destroy1(b);
}

TEST_CASE("get_end_motif homopolymer sequence") {
  bam1_t* b = make_read("AAAAAAAA", false);
  REQUIRE(feature_extractors::get_end_motif(b, 4) == "AAAA");
  bam_destroy1(b);
}
