import argparse
import csv
from pathlib import Path

import cfextract

from cfanalysis.classify import build_feature_matrix, run_loo_cv
from cfanalysis.features import metrics_to_features
from cfanalysis.plots import plot_end_motifs
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

    # --- Feature extraction ---
    samples: list[Sample] = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            bam_path = Path(row["bam_path"])
            if not bam_path.is_absolute():
                bam_path = manifest_dir / bam_path

            print(f"  {row['sample_id']} ({row['label']}) ...", end=" ", flush=True)
            metrics = cfextract.extract_features(str(bam_path), args.contig)
            features = metrics_to_features(metrics)
            samples.append(
                {
                    "sample_id": row["sample_id"],
                    "label": row["label"],
                    "features": features,
                }
            )
            print("done")

    print(f"\nProcessed {len(samples)} samples.")

    # --- Classification ---
    X, y, _ = build_feature_matrix(samples, top_n=args.top_motifs)
    result = run_loo_cv(X, y)

    print("\n--- LOO-CV Classification Report ---")
    print(result["report"])

    # --- Plot ---
    fig_path = output_dir / "end_motifs.png"
    plot_end_motifs(samples, top_n=10, out_path=fig_path)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
