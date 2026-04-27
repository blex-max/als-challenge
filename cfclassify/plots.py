from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.figure import Figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cfclassify.types import Sample


def plot_end_motifs(
    samples: list[Sample],
    top_n: int = 10,
    out_path: str | Path | None = None,
) -> Figure:
    """Grouped bar chart of top end-motif frequencies, ALS vs CTRL.

    Error bars show SEM across samples within each group.
    """
    groups: dict[str, dict[str, list[float]]] = {}
    for s in samples:
        label = s["label"]
        em = s["features"].get("end_motifs", {})
        if label not in groups:
            groups[label] = {}
        for motif, freq in em.items():
            groups[label].setdefault(motif, []).append(freq)

    all_freqs: dict[str, list[float]] = {}
    for per_label in groups.values():
        for motif, freqs in per_label.items():
            all_freqs.setdefault(motif, []).extend(freqs)
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
    ax.set_ylabel("Mean normalised frequency")
    ax.set_title(f"Top {top_n} end motifs — ALS vs CTRL (chr21, bisulfite cfDNA)")
    ax.legend()
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=150)

    return fig


def plot_frag_lengths(
    samples: list[Sample],
    out_path: str | Path | None = None,
) -> Figure:
    """Mean ± SEM fragment length distribution per group, 0–600 bp."""
    groups: dict[str, list[list[int]]] = {}
    for s in samples:
        label = s["label"]
        hist = s["features"].get("frag_len_hist", [])
        if hist:
            groups.setdefault(label, []).append(hist)

    x = np.arange(601)
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, hists in sorted(groups.items()):
        arr = np.array([h[:601] for h in hists], dtype=float)
        totals = arr.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1
        arr = arr / totals
        mean = arr.mean(axis=0)
        sem = arr.std(axis=0) / max(len(hists) ** 0.5, 1)
        ax.plot(x, mean, label=label.upper(), linewidth=1.2)
        ax.fill_between(x, mean - sem, mean + sem, alpha=0.2)

    ax.set_xlim(0, 600)
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Mean normalised frequency")
    ax.set_title("Fragment length distribution — ALS vs CTRL")
    ax.legend()
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=150)

    return fig


def plot_methylation(
    samples: list[Sample],
    out_path: str | Path | None = None,
) -> Figure:
    """Bar chart of mean CpG methylation rate per group with SEM and sample points."""
    groups: dict[str, list[float]] = {}
    for s in samples:
        rate = s["features"].get("methylation_mean")
        if rate is not None:
            groups.setdefault(s["label"], []).append(rate)

    sorted_labels = sorted(groups)
    x = np.arange(len(sorted_labels))
    means = [float(np.mean(groups[lb])) for lb in sorted_labels]
    sems = [
        float(np.std(groups[lb]) / max(len(groups[lb]) ** 0.5, 1))
        for lb in sorted_labels
    ]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.bar(x, means, yerr=sems, capsize=5, alpha=0.7)
    for i, lb in enumerate(sorted_labels):
        jitter = np.random.default_rng(0).uniform(-0.1, 0.1, len(groups[lb]))
        ax.scatter(i + jitter, groups[lb], color="black", s=20, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels([lb.upper() for lb in sorted_labels])
    ax.set_ylabel("Mean CpG methylation rate")
    ax.set_title("CpG methylation — ALS vs CTRL")
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=150)

    return fig


def plot_position_dist(
    samples: list[Sample],
    feature_key: str,
    title: str,
    out_path: str | Path | None = None,
) -> Figure:
    """Line plot of mean binned read position counts per group.

    feature_key should be 'start_pos_hist' or 'end_pos_hist'.
    Bin index represents position / 100_000 bp; x-axis shown in Mbp.
    """
    groups: dict[str, list[dict[int, int]]] = {}
    for s in samples:
        hist = s["features"].get(feature_key, {})
        if hist:
            groups.setdefault(s["label"], []).append(hist)

    if not groups:
        fig, ax = plt.subplots()
        ax.set_title(title)
        return fig

    all_bins: set[int] = set()
    for hists in groups.values():
        for h in hists:
            all_bins.update(h.keys())
    bins = sorted(all_bins)
    x_mbp = np.array([b * 0.1 for b in bins])

    fig, ax = plt.subplots(figsize=(12, 4))
    for label, hists in sorted(groups.items()):
        arr = np.array([[h.get(b, 0) for b in bins] for h in hists], dtype=float)
        mean = arr.mean(axis=0)
        sem = arr.std(axis=0) / max(len(hists) ** 0.5, 1)
        ax.plot(x_mbp, mean, label=label.upper(), linewidth=1.2)
        ax.fill_between(x_mbp, mean - sem, mean + sem, alpha=0.2)

    ax.set_xlabel("Position (Mbp)")
    ax.set_ylabel("Mean read count per 100 kbp bin")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=150)

    return fig
