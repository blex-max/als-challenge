#include "core/extract.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE (cfextract, m) {
  m.doc() = "Python bindings for cell-free DNA feature extractor core library";

  py::class_<cfextract::MethylationStats> (m, "MethylationStats")
    .def_readonly ("mean",      &cfextract::MethylationStats::mean)
    .def_readonly ("entropy",   &cfextract::MethylationStats::entropy)
    .def_readonly ("frac_high", &cfextract::MethylationStats::frac_high)
    .def_readonly ("frac_low",  &cfextract::MethylationStats::frac_low)
    .def_readonly ("median",    &cfextract::MethylationStats::median);

  py::class_<cfextract::FragLenStats> (m, "FragLenStats")
    .def_readonly ("mean",                &cfextract::FragLenStats::mean)
    .def_readonly ("std_dev",             &cfextract::FragLenStats::std_dev)
    .def_readonly ("frac_subnucleosomal", &cfextract::FragLenStats::frac_subnucleosomal)
    .def_readonly ("frac_nucleosomal",    &cfextract::FragLenStats::frac_nucleosomal)
    .def_readonly ("ratio_mono_di",       &cfextract::FragLenStats::ratio_mono_di);

  py::class_<cfextract::RegionMetrics, py::smart_holder> (m, "RegionMetrics")
    .def_readonly ("end_motifs",   &cfextract::RegionMetrics::end_motifs)
    .def_readonly ("frag_len_hist",&cfextract::RegionMetrics::frag_len_hist)
    .def_readonly ("fl_stats",     &cfextract::RegionMetrics::fl_stats)
    .def_readonly ("methylation",  &cfextract::RegionMetrics::methylation)
    .def_readonly ("start_pos_hist",&cfextract::RegionMetrics::start_pos_hist)
    .def_readonly ("end_pos_hist", &cfextract::RegionMetrics::end_pos_hist);

  // NOTE: coordinates are 0-based, half-open
  m.def ("extract_features",
    [](std::string fp, std::string contig, int64_t start, int64_t stop) {
      return cfextract::extract_metrics (fp, contig, start, stop);
    },
    py::arg ("aln_filepath"), py::arg ("contig"), py::arg ("region_start") = 0, py::arg ("region_end") = -1,
    "Extract cfDNA features from a BAM region. Contig is resolved by name from the BAM header."
  );
}
