#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace cfextract {

MethylationStats
compute_meth_stats (const std::unordered_map<int64_t, CpGSiteCount>& sites)
{
  std::vector<double> rates;
  rates.reserve (sites.size());
  size_t total_meth  = 0;
  size_t total_reads = 0;

  for (const auto& site : sites) {
    const size_t cov = site.second.methylated + site.second.unmethylated;
    if (cov == 0) { continue; }
    rates.push_back (static_cast<double> (site.second.methylated) / static_cast<double> (cov));
    total_meth  += site.second.methylated;
    total_reads += cov;
  }

  if (rates.empty()) { return {}; }

  const size_t n    = rates.size();
  const double mean = static_cast<double> (total_meth) / static_cast<double> (total_reads);

  double entropy = 0.0;
  size_t hi      = 0;
  size_t lo      = 0;
  for (const double r : rates) {
    // Binary entropy H(r) = -r log2(r) - (1-r) log2(1-r); 0 at boundaries by convention.
    if (r > 0.0 && r < 1.0) {
      entropy += -r * std::log2 (r) - (1.0 - r) * std::log2 (1.0 - r);
    }
    if (r > 0.8) { ++hi; }
    if (r < 0.1) { ++lo; }
  }
  entropy /= static_cast<double> (n);

  std::sort (rates.begin(), rates.end());

  return { mean, entropy,
           static_cast<double> (hi) / static_cast<double> (n),
           static_cast<double> (lo) / static_cast<double> (n),
           rates[n / 2] };
}


FragLenStats
compute_fl_stats (const std::array<size_t, 1001>& hist)
{
  const size_t total = std::accumulate (hist.begin(), hist.end(), size_t {0});
  if (total == 0) { return {}; }

  double mean = 0.0;
  for (size_t i = 0; i < hist.size(); ++i) {
    mean += static_cast<double> (i) * static_cast<double> (hist[i]);
  }
  mean /= static_cast<double> (total);

  double var = 0.0;
  for (size_t i = 0; i < hist.size(); ++i) {
    const double d = static_cast<double> (i) - mean;
    var += static_cast<double> (hist[i]) * d * d;
  }
  var /= static_cast<double> (total);

  // nucleosomal length boundaries follow the standard cfDNA literature
  const size_t sub  = std::accumulate (hist.begin(),       hist.begin() + 120, size_t {0});
  const size_t mono = std::accumulate (hist.begin() + 120, hist.begin() + 201, size_t {0});
  const size_t di   = std::accumulate (hist.begin() + 280, hist.begin() + 401, size_t {0});

  return { mean,
           std::sqrt (var),
           static_cast<double> (sub)  / static_cast<double> (total),
           static_cast<double> (mono) / static_cast<double> (total),
           static_cast<double> (mono) / static_cast<double> (std::max (di, size_t {1})) };
}

}
