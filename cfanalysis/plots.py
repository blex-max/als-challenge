from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.figure import Figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cfanalysis.types import Sample


def plot_end_motifs(
    samples: list[Sample],
    top_n: int = 10,
    out_path: str | Path | None = None,
) -> Figure:
    """Grouped bar chart of top end-motif frequencies, ALS vs CTRL.

    Error bars show SEM across samples within each group.
    """
    # Per-group frequency lists
    groups: dict[str, dict[str, list[float]]] = {}
    for s in samples:
        label = s["label"]
        em = s["features"].get("end_motifs", {})
        if label not in groups:
            groups[label] = {}
        for motif, freq in em.items():
            groups[label].setdefault(motif, []).append(freq)

    # Top-N motifs ranked by mean frequency across ALL samples
    all_freqs: dict[str, list[float]] = {}
    for s in samples:
        for motif, freq in s["features"].get("end_motifs", {}).items():
            all_freqs.setdefault(motif, []).append(freq)
    top_motifs = sorted(all_freqs, key=lambda m: -float(np.mean(all_freqs[m])))[:top_n]

    sorted_labels = sorted(groups)
    x = np.arange(len(top_motifs))
    bar_width = 0.8 / len(sorted_labels)

    fig, ax = plt.subplots(figsize=(12, 5))
    for i, label in enumerate(sorted_labels):
        means = [float(np.mean(groups[label].get(m, [0.0]))) for m in top_motifs]
        n = [len(groups[label].get(m, [0.0])) for m in top_motifs]
        sems = [
            float(np.std(groups[label].get(m, [0.0])) / max(n[j] ** 0.5, 1))
            for j, m in enumerate(top_motifs)
        ]
        offset = (i - len(sorted_labels) / 2 + 0.5) * bar_width
        ax.bar(x + offset, means, bar_width, label=label.upper(), yerr=sems, capsize=3)

    ax.set_xticks(x)
    ax.set_xticklabels(top_motifs, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("End motif (4-mer)")
    ax.set_ylabel("Mean normalized frequency")
    ax.set_title(f"Top {top_n} end motifs — ALS vs CTRL (chr21, bisulfite cfDNA)")
    ax.legend()
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=150)

    return fig
