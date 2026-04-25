from __future__ import annotations

from typing import TypedDict


class Features(TypedDict, total=False):
    end_motifs: dict[str, float]


class Sample(TypedDict):
    sample_id: str
    label: str
    features: Features
