import math

import cfextract

from cfanalysis.types import Features


def metrics_to_features(metrics: cfextract.RegionMetrics) -> Features:
    """Convert a RegionMetrics object into a plain Python feature dict.

    This is the sole extension point as RegionMetrics grows.
    Adding a new field means adding entries here; the rest of the pipeline
    picks up flat numeric features automatically.
    """
    features: Features = {}

    # End motifs: raw counts → normalised frequencies
    raw = dict(metrics.end_motifs)
    total_motifs = sum(raw.values())
    features["end_motifs"] = (
        {k: v / total_motifs for k, v in raw.items()} if total_motifs > 0 else {}
    )

    # Fragment length histogram (stored for plotting) and C++-computed summary stats
    features["frag_len_hist"] = list(metrics.frag_len_hist)
    fs = metrics.fl_stats
    if not math.isnan(fs.mean):
        features["fl_mean"] = fs.mean
        features["fl_std"] = fs.std_dev
        features["fl_frac_subnucleosomal"] = fs.frac_subnucleosomal
        features["fl_frac_nucleosomal"] = fs.frac_nucleosomal
        features["fl_ratio_mono_di"] = fs.ratio_mono_di

    # CpG methylation summary stats computed in C++ (NaN → absent keys when no coverage)
    ms = metrics.methylation
    if not math.isnan(ms.mean):
        features["methylation_mean"] = ms.mean
        features["methylation_entropy"] = ms.entropy
        features["methylation_frac_high"] = ms.frac_high
        features["methylation_frac_low"] = ms.frac_low
        features["methylation_median"] = ms.median

    # Position histograms (100 kbp bins) — for plotting only, not used in classification
    features["start_pos_hist"] = dict(metrics.start_pos_hist)
    features["end_pos_hist"] = dict(metrics.end_pos_hist)

    return features
