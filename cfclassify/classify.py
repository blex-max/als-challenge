import itertools

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import GridSearchCV, LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from cfclassify.types import Features, Sample


def generate_vocabulary(k: int) -> list[str]:
    return ["".join(b) for b in itertools.product("ACGT", repeat=k)]


def _make_inner_cv(n_train: int) -> StratifiedKFold:
    n_splits = max(2, min(5, n_train // 2))
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)


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


def run_loo_cv(samples: list[Sample], motif_k: int = 4) -> dict[str, object]:
    """Leave-one-out cross-validation with L2 logistic regression.

    Feature scaling and regularisation strength (C) are determined inside each
    fold on training samples only, preventing leakage from the held-out sample.
    C is chosen via inner StratifiedKFold GridSearchCV over [0.01, 0.1, 1.0, 10.0].
    The end-motif vocabulary is the full theoretical 4^k set — no per-fold selection.
    """
    predictions: list[str] = []
    all_labels = [s["label"] for s in samples]
    motifs = generate_vocabulary(motif_k)

    for train_idx, test_idx in LeaveOneOut().split(range(len(samples))):
        train_samples = [samples[i] for i in train_idx]
        test_samples = [samples[i] for i in test_idx]

        X_train, y_train, _ = build_feature_matrix(train_samples, motifs)
        X_test, _, _ = build_feature_matrix(test_samples, motifs)

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        np.nan_to_num(X_train_s, copy=False)
        np.nan_to_num(X_test_s, copy=False)

        grid = GridSearchCV(
            LogisticRegression(max_iter=1000, random_state=42),
            {"C": [0.01, 0.1, 1.0, 10.0]},
            cv=_make_inner_cv(len(train_samples)),
            refit=True,
        )
        grid.fit(X_train_s, y_train)
        predictions.append(str(grid.predict(X_test_s)[0]))

    return {
        "predictions": predictions,
        "labels": all_labels,
        "report": classification_report(all_labels, predictions),
    }


def train_final_model(
    samples: list[Sample],
    motif_k: int = 4,
) -> tuple[LogisticRegression, StandardScaler, list[str]]:
    """Train a final L2 logistic regression on ALL samples.

    C is selected by inner stratified k-fold CV. Returns (model, scaler, flat_keys);
    flat_keys preserves the non-motif feature ordering needed for inference alignment.
    """
    motifs = generate_vocabulary(motif_k)
    X, y, feature_names = build_feature_matrix(samples, motifs)
    flat_keys = [fn for fn in feature_names if not fn.startswith("motif_")]
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    np.nan_to_num(X_s, copy=False)
    grid = GridSearchCV(
        LogisticRegression(max_iter=1000, random_state=42),
        {"C": [0.01, 0.1, 1.0, 10.0]},
        cv=_make_inner_cv(len(samples)),
        refit=True,
    )
    grid.fit(X_s, y)
    model: LogisticRegression = grid.best_estimator_
    return model, scaler, flat_keys


def predict_sample(
    sample_features: Features,
    model: LogisticRegression,
    scaler: StandardScaler,
    flat_keys: list[str],
    motif_k: int = 4,
) -> tuple[str, float]:
    """Predict class label and confidence for a single sample's feature dict."""
    motifs = generate_vocabulary(motif_k)
    em = sample_features.get("end_motifs", {})
    motif_vals = [em.get(m, 0.0) for m in motifs]
    flat_vals = [sample_features.get(k, float("nan")) for k in flat_keys]  # type: ignore[literal-required]
    x = np.array([motif_vals + flat_vals], dtype=float)
    x_s = scaler.transform(x)
    np.nan_to_num(x_s, copy=False)
    label: str = str(model.predict(x_s)[0])
    proba: float = float(model.predict_proba(x_s)[0].max())
    return label, proba
