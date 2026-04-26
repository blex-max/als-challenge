from typing import TypedDict


class Features(TypedDict, total=False):
    end_motifs: dict[str, float]
    frag_len_hist: list[int]
    fl_mean: float
    fl_std: float
    fl_frac_subnucleosomal: float
    fl_frac_nucleosomal: float
    fl_ratio_mono_di: float
    methylation_mean: float
    methylation_entropy: float
    methylation_frac_high: float
    methylation_frac_low: float
    methylation_median: float
    start_pos_hist: dict[int, int]
    end_pos_hist: dict[int, int]


class Sample(TypedDict):
    sample_id: str
    label: str
    features: Features
