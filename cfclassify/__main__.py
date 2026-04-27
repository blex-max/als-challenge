import argparse
import csv
import time
from collections.abc import Callable
from pathlib import Path

import cfextract
import matplotlib.pyplot as plt

from cfclassify.classify import (
    predict_sample,
    run_loo_cv,
    train_final_model,
)
from cfclassify.features import metrics_to_features
from cfclassify.model import (
    load_feature_cache,
    load_model,
    save_feature_cache,
    save_model,
)
from cfclassify.plots import (
    plot_end_motifs,
    plot_frag_lengths,
    plot_methylation,
    plot_position_dist,
)
from cfclassify.types import ModelBundle, Sample


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cfclassify",
        description=(
            "cfDNA analysis pipeline — end-motif extraction and ALS/CTRL classification"
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_train = sub.add_parser(
        "train",
        help="Train a model for classifying samples on a labelled manifest",
    )
    p_train.add_argument(
        "--manifest",
        required=True,
        help="CSV file with columns: sample_id, bam_path, label",
    )
    p_train.add_argument(
        "--contig", default="chr21", help="Contig name to process (default: chr21)"
    )
    p_train.add_argument(
        "--out-dir",
        default=".",
        help=(
            "Directory for saving trained model and model evaluation"
            " results (default: ./)"
        ),
    )
    p_train.add_argument(
        "--motif-length",
        type=int,
        default=4,
        help="k-mer length for end-motif features; 4^k features used (default: 4)",
    )

    p_predict = sub.add_parser(
        "predict", help="Predict class for a single unlabelled BAM"
    )
    p_predict.add_argument("--bam", required=True, help="Path to the BAM file")
    p_predict.add_argument(
        "--model-path", required=True, help="Path to a saved model bundle"
    )
    p_predict.add_argument(
        "--sample-id",
        default="unknown",
        help="Label for progress output (default: unknown)",
    )

    p_update = sub.add_parser(
        "update",
        help="Add a new labelled sample to the feature cache and retrain the model",
    )
    p_update.add_argument("--bam", required=True, help="BAM file for the new sample")
    p_update.add_argument(
        "--label", required=True, help="Class label for the new sample"
    )
    p_update.add_argument(
        "--sample-id", required=True, help="Unique identifier for the new sample"
    )
    p_update.add_argument(
        "--model-path", required=True, help="Path to an existing model bundle"
    )

    args = parser.parse_args()
    dispatch: dict[str, Callable[[argparse.Namespace], None]] = {
        "train": _run_train,
        "predict": _run_predict,
        "update": _run_update,
    }
    dispatch[args.command](args)


def _load_samples(
    manifest_path: Path,
    contig: str,
    motif_k: int,
) -> list[Sample]:
    manifest_dir = manifest_path.parent
    samples: list[Sample] = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            bam_path = Path(row["bam_path"])
            if not bam_path.is_absolute():
                bam_path = manifest_dir / bam_path
            print(f"  {row['sample_id']} ({row['label']}) ...", end=" ", flush=True)
            metrics = cfextract.extract_features(
                str(bam_path), contig, motif_sz=motif_k
            )
            features = metrics_to_features(metrics)
            print("done")
            samples.append(
                {
                    "sample_id": row["sample_id"],
                    "label": row["label"],
                    "features": features,
                }
            )
    return samples


def _run_train(args: argparse.Namespace) -> None:
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = _load_samples(Path(args.manifest), args.contig, args.motif_length)
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

    result = run_loo_cv(samples, motif_k=args.motif_length)

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

    model, scaler, flat_keys = train_final_model(samples, motif_k=args.motif_length)
    bundle: ModelBundle = {
        "model": model,
        "scaler": scaler,
        "k": args.motif_length,
        "flat_keys": flat_keys,
        "contig": args.contig,
        "n_training_samples": len(samples),
    }
    model_path = output_dir / "model.pkl"
    save_model(bundle, model_path)
    cache_path = save_feature_cache(samples, model_path)
    print(f"Saved model:         {model_path}")
    print(f"Saved feature cache: {cache_path}")


def _run_predict(args: argparse.Namespace) -> None:
    bundle = load_model(Path(args.model_path))
    print(
        f"Loaded model trained on {bundle['n_training_samples']} samples "
        f"(k={bundle['k']}, contig={bundle['contig']})"
    )
    print(f"  {args.sample_id} ...", end=" ", flush=True)
    t0 = time.perf_counter()
    metrics = cfextract.extract_features(
        args.bam, bundle["contig"], motif_sz=bundle["k"]
    )
    elapsed = time.perf_counter() - t0
    features = metrics_to_features(metrics)
    print(f"done ({elapsed:.1f}s)")
    label, proba = predict_sample(
        features,
        bundle["model"],
        bundle["scaler"],
        bundle["flat_keys"],
        motif_k=bundle["k"],
    )
    print(f"Prediction: {label}  (p={proba:.3f})")


def _run_update(args: argparse.Namespace) -> None:
    model_path = Path(args.model_path)
    bundle = load_model(model_path)

    print(f"  {args.sample_id} ({args.label}) ...", end=" ", flush=True)
    t0 = time.perf_counter()
    metrics = cfextract.extract_features(
        args.bam, bundle["contig"], motif_sz=bundle["k"]
    )
    elapsed = time.perf_counter() - t0
    features = metrics_to_features(metrics)
    print(f"done ({elapsed:.1f}s)")

    new_sample: Sample = {
        "sample_id": args.sample_id,
        "label": args.label,
        "features": features,
    }

    samples = load_feature_cache(model_path)
    samples.append(new_sample)
    _validate_cache(samples)

    model, scaler, flat_keys = train_final_model(samples, motif_k=bundle["k"])
    new_bundle: ModelBundle = {
        "model": model,
        "scaler": scaler,
        "k": bundle["k"],
        "flat_keys": flat_keys,
        "contig": bundle["contig"],
        "n_training_samples": len(samples),
    }
    save_model(new_bundle, model_path)
    save_feature_cache(samples, model_path)
    print(f"Updated model saved to {model_path} ({len(samples)} total samples)")


def _validate_cache(samples: list[Sample]) -> None:
    """Raise ValueError if cache is too small for meaningful retraining."""
    labels = [s["label"] for s in samples]
    classes = set(labels)
    if len(classes) < 2:
        raise ValueError(
            f"Cache contains only one class ({classes!r}). "
            "Add at least one sample from each class before retraining."
        )
    min_count = min(labels.count(c) for c in classes)
    if min_count < 2:
        raise ValueError(
            f"Each class needs at least 2 samples for inner CV. "
            f"Smallest class has {min_count} sample(s)."
        )


if __name__ == "__main__":
    main()
