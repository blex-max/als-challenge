#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "core/core.hpp"

namespace py = pybind11;
namespace cfe = cfextract;

PYBIND11_MODULE (cfextract, m) {
  m.doc() = "Python bindings for cell-free DNA feature extractor core library";

  py::class_<cfe::MethylationStats> (m, "MethylationStats")
    .def_readonly ("mean", &cfe::MethylationStats::mean)
    .def_readonly ("entropy", &cfe::MethylationStats::entropy)
    .def_readonly ("frac_high", &cfe::MethylationStats::frac_high)
    .def_readonly ("frac_low", &cfe::MethylationStats::frac_low)
    .def_readonly ("median", &cfe::MethylationStats::median);

  py::class_<cfe::FragLenStats> (m, "FragLenStats")
    .def_readonly ("mean", &cfe::FragLenStats::mean)
    .def_readonly ("std_dev", &cfe::FragLenStats::std_dev)
    .def_readonly ("frac_subnucleosomal", &cfe::FragLenStats::frac_subnucleosomal)
    .def_readonly ("frac_nucleosomal", &cfe::FragLenStats::frac_nucleosomal)
    .def_readonly ("ratio_mono_di", &cfe::FragLenStats::ratio_mono_di);

  py::class_<cfe::RegionMetrics, py::smart_holder> (m, "RegionMetrics")
    .def_readonly ("end_motifs", &cfe::RegionMetrics::end_motifs)
    .def_readonly ("frag_len_hist", &cfe::RegionMetrics::frag_len_hist)
    .def_readonly ("fl_stats", &cfe::RegionMetrics::fl_stats)
    .def_readonly ("methylation", &cfe::RegionMetrics::methylation)
    .def_readonly ("start_pos_hist", &cfe::RegionMetrics::start_pos_hist)
    .def_readonly ("end_pos_hist", &cfe::RegionMetrics::end_pos_hist);

  // NOTE: coordinates are 0-based, half-open
  m.def ("extract_features",
    [](std::string fp, std::string contig, int64_t start, int64_t stop, size_t motif_sz) {
      return cfe::extract_metrics (fp, contig, start, stop, motif_sz);
    },
    py::arg ("aln_filepath"),
    py::arg ("contig"),
    py::arg ("region_start") = 0,
    py::arg ("region_end") = -1,
    py::arg ("motif_sz") = size_t{4},
    "Extract cfDNA features from a BAM region. Contig is resolved by name from the BAM header."
  );
}
