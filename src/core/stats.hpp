#pragma once

#include <array>
#include <cstddef>
#include <limits>

#include "extract.hpp"

namespace cfextract {

struct MethylationStats {
  double mean = std::numeric_limits<double>::quiet_NaN();
  double entropy = std::numeric_limits<double>::quiet_NaN();
  double frac_high = std::numeric_limits<double>::quiet_NaN();
  double frac_low = std::numeric_limits<double>::quiet_NaN();
  double median = std::numeric_limits<double>::quiet_NaN();
};

struct FragLenStats {
  double mean = std::numeric_limits<double>::quiet_NaN();
  double std_dev = std::numeric_limits<double>::quiet_NaN();
  double frac_subnucleosomal = std::numeric_limits<double>::quiet_NaN();
  double frac_nucleosomal = std::numeric_limits<double>::quiet_NaN();
  double ratio_mono_di = std::numeric_limits<double>::quiet_NaN();
};


MethylationStats
compute_meth_stats
(const extractors::CpGMap& sites);


// fragment lengths: index = length in bp, value = count; >1000 bp clamped to index 1000
using FragLenBins = std::array<size_t, 1001>;

FragLenStats
compute_fl_stats
(const FragLenBins& hist);

}
