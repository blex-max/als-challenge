#pragma once

#include <string>

#include <htslib/sam.h>

// individual feature extractors
namespace feature_extractors {

std::string
get_end_motif
(const bam1_t* b, size_t motif_sz);

}
