import json
from pathlib import Path
from typing import cast

import joblib
import numpy as np

from cfclassify.types import Features, ModelBundle, Sample


def save_model(bundle: ModelBundle, path: Path) -> None:
    joblib.dump(bundle, path)


def load_model(path: Path) -> ModelBundle:
    bundle: ModelBundle = joblib.load(path)
    return bundle


def save_feature_cache(samples: list[Sample], model_path: Path) -> Path:
    """Write feature cache as JSON alongside model_path; return the cache path."""
    cache_path = _cache_path(model_path)
    records = [
        {
            "sample_id": s["sample_id"],
            "label": s["label"],
            "features": _serialise_features(s["features"]),
        }
        for s in samples
    ]
    cache_path.write_text(json.dumps(records, indent=2))
    return cache_path


def load_feature_cache(model_path: Path) -> list[Sample]:
    """Load feature cache JSON from alongside model_path."""
    cache_path = _cache_path(model_path)
    records: list[dict[str, object]] = json.loads(cache_path.read_text())
    return [
        {
            "sample_id": str(r["sample_id"]),
            "label": str(r["label"]),
            "features": _deserialise_features(cast(dict[str, object], r["features"])),
        }
        for r in records
    ]


def _cache_path(model_path: Path) -> Path:
    return model_path.with_suffix(model_path.suffix + ".features.json")


def _serialise_features(features: Features) -> dict[str, object]:
    out: dict[str, object] = {}
    for k, v in features.items():
        if isinstance(v, np.integer):
            out[k] = int(v)
        else:
            out[k] = v
    return out


def _deserialise_features(raw: dict[str, object]) -> Features:
    features: Features = {}
    for k, v in raw.items():
        if k == "end_motifs":
            d = cast(dict[str, float], v)
            features["end_motifs"] = {str(mk): float(mv) for mk, mv in d.items()}
        elif k == "frag_len_hist":
            features["frag_len_hist"] = cast(list[int], v)
        elif k in ("start_pos_hist", "end_pos_hist"):
            d2 = cast(dict[str, int], v)
            features[k] = {int(bk): int(bv) for bk, bv in d2.items()}  # type: ignore[literal-required]
        elif isinstance(v, (int, float)):
            features[k] = float(v)  # type: ignore[literal-required]
    return features
