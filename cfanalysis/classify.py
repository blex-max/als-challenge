from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from cfanalysis.types import Sample


def build_feature_matrix(
    samples: list[Sample],
    top_n: int = 20,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Build (X, y, feature_names) from a list of sample dicts.

    Each dict has keys: sample_id, label, features.
    Top-N end motifs are selected by summed frequency across all samples.
    Any flat numeric keys in features (e.g. methylation_rate, fl_mean added
    in future) are appended automatically.
    """
    # Top-N end motifs by total frequency across all samples
    motif_totals: dict[str, float] = {}
    for s in samples:
        for motif, freq in s["features"].get("end_motifs", {}).items():
            motif_totals[motif] = motif_totals.get(motif, 0.0) + freq
    top_motifs = sorted(motif_totals, key=lambda m: -motif_totals[m])[:top_n]

    # Any flat numeric features present alongside end_motifs
    flat_keys = (
        [
            k
            for k in samples[0]["features"]
            if k != "end_motifs" and isinstance(samples[0]["features"][k], (int, float))
        ]
        if samples
        else []
    )

    feature_names = [f"motif_{m}" for m in top_motifs] + flat_keys

    rows = []
    labels = []
    for s in samples:
        em = s["features"].get("end_motifs", {})
        motif_vals = [em.get(m, 0.0) for m in top_motifs]
        flat_vals = [s["features"].get(k, float("nan")) for k in flat_keys]
        rows.append(motif_vals + flat_vals)
        labels.append(s["label"])

    return np.array(rows, dtype=float), np.array(labels), feature_names


def run_loo_cv(X: np.ndarray, y: np.ndarray) -> dict[str, object]:
    """Leave-one-out cross-validation with L2 logistic regression.

    StandardScaler and the classifier are both fit inside each fold to
    avoid data leakage.
    """
    predictions: list[str] = []

    for train_idx, test_idx in LeaveOneOut().split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        np.nan_to_num(X_train_s, copy=False)
        np.nan_to_num(X_test_s, copy=False)

        clf = LogisticRegression(C=0.1, max_iter=1000, random_state=42)
        clf.fit(X_train_s, y_train)
        predictions.append(str(clf.predict(X_test_s)[0]))

    return {
        "predictions": predictions,
        "labels": y.tolist(),
        "report": classification_report(y, predictions),
    }
