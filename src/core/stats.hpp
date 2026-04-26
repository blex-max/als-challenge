#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

#include "extract.hpp"

namespace cfextract {

MethylationStats
compute_meth_stats (const std::unordered_map<int64_t, CpGSiteCount>& sites);

FragLenStats
compute_fl_stats (const std::array<size_t, 1001>& hist);

}
