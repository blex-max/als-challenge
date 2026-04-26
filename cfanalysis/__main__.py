import argparse
import csv
import time
from pathlib import Path

import cfextract
import matplotlib.pyplot as plt

from cfanalysis.classify import run_loo_cv
from cfanalysis.features import metrics_to_features
from cfanalysis.plots import (
    plot_end_motifs,
    plot_frag_lengths,
    plot_methylation,
    plot_position_dist,
)
from cfanalysis.types import Sample


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cfanalysis",
        description=(
            "cfDNA analysis pipeline — end-motif extraction and ALS/CTRL classification"
        ),
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="CSV file with columns: sample_id, bam_path, label",
    )
    parser.add_argument(
        "--contig", default="chr21", help="Contig name to process (default: chr21)"
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory for output files (default: .)"
    )
    parser.add_argument(
        "--top-motifs",
        type=int,
        default=20,
        help="Number of top end motifs to use as classification features (default: 20)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)
    manifest_dir = manifest_path.parent

    samples: list[Sample] = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            bam_path = Path(row["bam_path"])
            if not bam_path.is_absolute():
                bam_path = manifest_dir / bam_path

            print(f"  {row['sample_id']} ({row['label']}) ...", end=" ", flush=True)
            t0 = time.perf_counter()
            metrics = cfextract.extract_features(str(bam_path), args.contig)
            elapsed = time.perf_counter() - t0
            features = metrics_to_features(metrics)
            samples.append(
                {
                    "sample_id": row["sample_id"],
                    "label": row["label"],
                    "features": features,
                }
            )
            print(f"done ({elapsed:.1f}s)")

    print(f"\nProcessed {len(samples)} samples.")

    summary_fields = [
        "fl_mean",
        "fl_std",
        "fl_frac_subnucleosomal",
        "fl_frac_nucleosomal",
        "fl_ratio_mono_di",
        "methylation_mean",
        "methylation_entropy",
        "methylation_frac_high",
        "methylation_frac_low",
        "methylation_median",
    ]
    summary_path = output_dir / "sample_summary.csv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample_id", "label"] + summary_fields, extrasaction="ignore"
        )
        writer.writeheader()
        for s in samples:
            out_row: dict[str, object] = {
                "sample_id": s["sample_id"],
                "label": s["label"],
            }
            for field in summary_fields:
                out_row[field] = s["features"].get(field, "")  # type: ignore[literal-required]
            writer.writerow(out_row)
    print(f"Saved: {summary_path}")

    result = run_loo_cv(samples, top_n=args.top_motifs)

    report_str = str(result["report"])
    report_path = output_dir / "classification_report.txt"
    report_path.write_text(report_str)
    print(f"Saved: {report_path}")

    print("\n--- LOO-CV Classification Report ---")
    print(report_str)

    for name, fig in [
        ("end_motifs.png", plot_end_motifs(samples, top_n=10)),
        ("frag_lengths.png", plot_frag_lengths(samples)),
        ("methylation.png", plot_methylation(samples)),
        (
            "start_positions.png",
            plot_position_dist(
                samples, "start_pos_hist", "Fragment start position distribution"
            ),
        ),
        (
            "end_positions.png",
            plot_position_dist(
                samples, "end_pos_hist", "Fragment end position distribution"
            ),
        ),
    ]:
        path = output_dir / name
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
