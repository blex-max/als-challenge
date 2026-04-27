#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include <htslib/sam.h>

#include "core/access.hpp"
#include "core/extract.hpp"
#include "core/stats.hpp"

// minimal synthetic bam1_t with known sequence
static bam1_t*
make_read (const std::string& seq, bool reverse_strand)
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
    if (i % 2 == 0) {
      packed[i / 2]  = static_cast<uint8_t>(code << 4U);
    }
    else {
      packed[i / 2] |= code;
    }
  }
  return b;
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
// accumlate_motifs
// ---------------------------------------------------------------------------

TEST_CASE("accumlate_motifs forward strand records first N bases") {
  bam1_t* b = make_read("ACGTTTTT", false);
  extractors::MotifCounts counts;
  extractors::accumlate_motifs(counts, b, 4);
  REQUIRE(counts["ACGT"] == 1);
  bam_destroy1(b);
}

TEST_CASE("accumlate_motifs reverse strand records reverse complement of last N bases") {
  // seq: TTTTACGA  last 4: ACGA  reversed: AGCA  complemented: TCGT
  bam1_t* b = make_read("TTTTACGA", true);
  extractors::MotifCounts counts;
  extractors::accumlate_motifs(counts, b, 4);
  REQUIRE(counts["TCGT"] == 1);
  bam_destroy1(b);
}

TEST_CASE("accumlate_motifs respects smaller motif size") {
  bam1_t* b = make_read("ACGTTTTT", false);
  extractors::MotifCounts counts;
  extractors::accumlate_motifs(counts, b, 2);
  REQUIRE(counts["AC"] == 1);
  bam_destroy1(b);
}

TEST_CASE("accumlate_motifs homopolymer sequence") {
  bam1_t* b = make_read("AAAAAAAA", false);
  extractors::MotifCounts counts;
  extractors::accumlate_motifs(counts, b, 4);
  REQUIRE(counts["AAAA"] == 1);
  bam_destroy1(b);
}

TEST_CASE("accumlate_motifs reverse strand homopolymer is self-complementary") {
  // RC of TTTT is AAAA
  bam1_t* b = make_read("AAAATTTT", true);
  extractors::MotifCounts counts;
  extractors::accumlate_motifs(counts, b, 4);
  REQUIRE(counts["AAAA"] == 1);
  bam_destroy1(b);
}

// ---------------------------------------------------------------------------
// compute_fl_stats
// ---------------------------------------------------------------------------

TEST_CASE("compute_fl_stats empty histogram returns NaN") {
  std::array<size_t, 1001> hist{};
  const auto s = cfextract::compute_fl_stats (hist);
  REQUIRE(std::isnan(s.mean));
  REQUIRE(std::isnan(s.std_dev));
}

TEST_CASE("compute_fl_stats all fragments at one length: mean correct, std zero") {
  std::array<size_t, 1001> hist{};
  hist[150] = 1000;
  const auto s = cfextract::compute_fl_stats (hist);
  REQUIRE(s.mean == Catch::Approx(150.0));
  REQUIRE(s.std_dev == Catch::Approx(0.0).margin(1e-9));
  // 150 bp falls in the mono-nucleosomal window (120–200), so all fragments are nucleosomal
  REQUIRE(s.frac_nucleosomal == Catch::Approx(1.0));
  REQUIRE(s.frac_subnucleosomal == Catch::Approx(0.0));
}

TEST_CASE("compute_fl_stats mono/di ratio is 1 when counts are equal") {
  std::array<size_t, 1001> hist{};
  hist[150] = 500;   // mono-nucleosomal bin
  hist[320] = 500;   // di-nucleosomal bin
  const auto s = cfextract::compute_fl_stats (hist);
  REQUIRE(s.ratio_mono_di == Catch::Approx(1.0));
}

// ---------------------------------------------------------------------------
// compute_meth_stats
// ---------------------------------------------------------------------------

TEST_CASE("compute_meth_stats empty site map returns NaN") {
  extractors::CpGMap sites;
  const auto s = cfextract::compute_meth_stats (sites);
  REQUIRE(std::isnan(s.mean));
  REQUIRE(std::isnan(s.entropy));
}

TEST_CASE("compute_meth_stats all-methylated sites: mean 1, entropy 0, frac_high 1") {
  extractors::CpGMap sites;
  for (int64_t i = 0; i < 10; ++i) {
    sites[i] = { 10, 0 };   // methylated=10, unmethylated=0 → rate=1.0
  }
  const auto s = cfextract::compute_meth_stats (sites);
  REQUIRE(s.mean      == Catch::Approx(1.0));
  REQUIRE(s.entropy   == Catch::Approx(0.0).margin(1e-9));
  REQUIRE(s.frac_high == Catch::Approx(1.0));
  REQUIRE(s.frac_low  == Catch::Approx(0.0));
  REQUIRE(s.median    == Catch::Approx(1.0));
}

TEST_CASE("compute_meth_stats all-unmethylated sites: mean 0, frac_low 1") {
  extractors::CpGMap sites;
  for (int64_t i = 0; i < 10; ++i) {
    sites[i] = { 0, 10 };
  }
  const auto s = cfextract::compute_meth_stats (sites);
  REQUIRE(s.mean      == Catch::Approx(0.0).margin(1e-9));
  REQUIRE(s.entropy   == Catch::Approx(0.0).margin(1e-9));
  REQUIRE(s.frac_high == Catch::Approx(0.0).margin(1e-9));
  REQUIRE(s.frac_low  == Catch::Approx(1.0));
}
