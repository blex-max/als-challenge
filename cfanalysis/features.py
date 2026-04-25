from __future__ import annotations

import cfextract

from cfanalysis.types import Features


def metrics_to_features(metrics: cfextract.RegionMetrics) -> Features:
    """Convert a RegionMetrics object into a plain Python feature dict.

    This is the sole extension point as RegionMetrics grows.
    Adding a new field means adding a new key here; the rest of the
    pipeline picks it up automatically.
    """
    features: Features = {}

    # End motifs: raw counts → normalized frequencies
    raw = dict(metrics.end_motifs)
    total = sum(raw.values())
    features["end_motifs"] = {k: v / total for k, v in raw.items()} if total > 0 else {}

    # Future fields (uncomment as RegionMetrics gains them):
    # features["fl_mean"]          = ...
    # features["fl_median"]        = ...
    # features["fl_frac_nfr"]      = ...
    # features["methylation_rate"] = ...

    return features
