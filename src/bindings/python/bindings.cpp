#include "core/extract.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE (cfextract, m) {
  m.doc() = "Python bindings for cell-free DNA feature extractor core library";

  py::class_<cfextract::RegionMetrics, py::smart_holder> (m, "RegionMetrics")
    .def_readonly ("end_motifs", &cfextract::RegionMetrics::end_motifs);

  // NOTE: must document coordinate conventions (0-based, half-open)
  m.def ("extract_features",
    [](std::string fp, std::string contig, int64_t start, int64_t stop) {
      return cfextract::extract_metrics (fp, contig, start, stop);
    },
    py::arg ("aln_filepath"), py::arg ("contig"), py::arg ("region_start") = 0, py::arg ("region_end") = -1,
    "Extract cfDNA features from a BAM region. Contig is resolved by name from the BAM header."
  );
}
