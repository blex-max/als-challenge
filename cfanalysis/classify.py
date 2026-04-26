import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import GridSearchCV, LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from cfanalysis.types import Sample


def _select_top_motifs(samples: list[Sample], top_n: int) -> list[str]:
    totals: dict[str, float] = {}
    for s in samples:
        for motif, freq in s["features"].get("end_motifs", {}).items():
            totals[motif] = totals.get(motif, 0.0) + freq
    return sorted(totals, key=lambda m: -totals[m])[:top_n]


def build_feature_matrix(
    samples: list[Sample],
    motifs: list[str],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Build (X, y, feature_names) from a list of samples and an explicit motif list.

    Flat numeric keys in features (fl_mean, methylation_mean, etc.) are appended
    automatically. Non-numeric keys (dicts, lists) are excluded by the isinstance check.
    """
    flat_keys = (
        [
            k
            for k in samples[0]["features"]
            if k != "end_motifs" and isinstance(samples[0]["features"][k], (int, float))
        ]
        if samples
        else []
    )

    feature_names = [f"motif_{m}" for m in motifs] + flat_keys

    rows = []
    labels = []
    for s in samples:
        em = s["features"].get("end_motifs", {})
        motif_vals = [em.get(m, 0.0) for m in motifs]
        flat_vals = [s["features"].get(k, float("nan")) for k in flat_keys]  # type: ignore[literal-required]
        rows.append(motif_vals + flat_vals)
        labels.append(s["label"])

    return np.array(rows, dtype=float), np.array(labels), feature_names


def run_loo_cv(samples: list[Sample], top_n: int = 20) -> dict[str, object]:
    """Leave-one-out cross-validation with L2 logistic regression.

    End-motif selection, feature scaling, and regularisation strength (C) are
    all determined inside each fold on training samples only, preventing leakage
    from the held-out sample. C is chosen via inner 5-fold StratifiedKFold
    GridSearchCV over [0.01, 0.1, 1.0, 10.0].
    """
    predictions: list[str] = []
    all_labels = [s["label"] for s in samples]

    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    param_grid = {"C": [0.01, 0.1, 1.0, 10.0]}

    for train_idx, test_idx in LeaveOneOut().split(range(len(samples))):
        train_samples = [samples[i] for i in train_idx]
        test_samples = [samples[i] for i in test_idx]

        motifs = _select_top_motifs(train_samples, top_n)
        X_train, y_train, _ = build_feature_matrix(train_samples, motifs)
        X_test, _, _ = build_feature_matrix(test_samples, motifs)

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        np.nan_to_num(X_train_s, copy=False)
        np.nan_to_num(X_test_s, copy=False)

        grid = GridSearchCV(
            LogisticRegression(max_iter=1000, random_state=42),
            param_grid,
            cv=inner_cv,
            refit=True,
        )
        grid.fit(X_train_s, y_train)
        predictions.append(str(grid.predict(X_test_s)[0]))

    return {
        "predictions": predictions,
        "labels": all_labels,
        "report": classification_report(all_labels, predictions),
    }
