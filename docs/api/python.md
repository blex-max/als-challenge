# Python API — cfclassify

All public symbols live in the `cfclassify` package. Import paths are shown on each entry.

---

## Types — `cfclassify.types`

Three `TypedDict` classes define the data schema flowing through the Python layer.

```python
class Features(TypedDict, total=False):
    end_motifs: dict[str, float]         # normalised k-mer frequencies
    frag_len_hist: list[int]             # 1001-bin histogram (0–1000 bp)
    fl_mean: float
    fl_std: float
    fl_frac_subnucleosomal: float        # fraction < 120 bp
    fl_frac_nucleosomal: float           # fraction 120–200 bp
    fl_ratio_mono_di: float              # mono / di-nucleosomal count ratio
    methylation_mean: float              # coverage-weighted mean
    methylation_entropy: float           # mean binary entropy across sites
    methylation_frac_high: float         # fraction of sites with rate > 0.8
    methylation_frac_low: float          # fraction of sites with rate < 0.1
    methylation_median: float
    start_pos_hist: dict[int, int]       # 100 kbp bins, for plotting only
    end_pos_hist: dict[int, int]

class Sample(TypedDict):
    sample_id: str
    label: str
    features: Features

class ModelBundle(TypedDict):
    model: LogisticRegression
    scaler: StandardScaler
    k: int                    # motif length used during training
    flat_keys: list[str]      # non-motif feature names, for vector alignment
    contig: str
    n_training_samples: int
```

**`Features`** is `total=False`: every key is optional. Missing data (e.g. no CpG coverage in a region) is represented by the absence of the corresponding key — not by a `NaN` value. This keeps downstream code honest: a consumer must explicitly handle the absent-key case rather than silently propagating NaN.

**`Sample`** wraps a single processed BAM file.

**`ModelBundle`** holds everything needed to classify a new sample: the fitted model, fitted scaler, k-mer length `k`, the ordered list of non-motif feature names `flat_keys` (for vector alignment at inference time), and the contig and sample count used during training.

---

## Feature conversion — `cfclassify.features`

```python
def metrics_to_features(metrics: cfextract.RegionMetrics) -> Features: ...
```

`metrics_to_features()` is the **sole extension point** for adding new features to the pipeline. When a new field is added to `RegionMetrics` on the C++ side:

1. Add the conversion logic here.
2. If the new field is a flat numeric scalar, `build_feature_matrix()` in `classify.py` picks it up automatically — no other changes needed.
3. If the new field is a histogram or dict (like `end_motifs`), add explicit handling to both `metrics_to_features()` and `build_feature_matrix()`.

---

## Classification — `cfclassify.classify`

```python
def generate_vocabulary(k: int) -> list[str]: ...
def run_loo_cv(samples: list[Sample], motif_k: int = 4) -> dict[str, object]: ...
def train_final_model(samples: list[Sample], motif_k: int = 4) -> tuple[...]: ...
def predict_sample(sample_features: Features, model, scaler, flat_keys, motif_k: int = 4) -> tuple[str, float]: ...
def build_feature_matrix(samples: list[Sample], motifs: list[str]) -> tuple[...]: ...
```

### `generate_vocabulary(k)`

Returns the full sorted k-mer vocabulary as a `list[str]` with `4^k` entries over `{A, C, G, T}`, generated via `itertools.product`. This is the vocabulary used as classifier features — deterministic from `k` alone, no training data required.

### `run_loo_cv(samples, motif_k=4)`

Runs leave-one-out cross-validation and returns a dict:

```python
{
    "predictions": list[str],   # one label per sample, in input order
    "labels":      list[str],   # ground-truth labels in the same order
    "report":      str,         # sklearn classification_report string
}
```

The two leakage guards (scaling, hyperparameter search) are applied inside the fold loop. See [Classification — Implementation](../classification.md#implementation) for an annotated walkthrough.

### `train_final_model(samples, motif_k=4)`

Trains a final L2 logistic regression on all samples (not LOO-CV). Returns `(model, scaler, flat_keys)` — the three components needed to build a `ModelBundle` for deployment.

### `predict_sample(sample_features, model, scaler, flat_keys, motif_k=4)`

Classifies a single `Features` dict using a pre-trained model bundle. Returns `(label, probability)` where `probability` is the confidence for the predicted class.

### `build_feature_matrix(samples, motifs)`

Builds `(X, y, feature_names)` from a sample list and an explicit motif list. Flat numeric keys in `Features` (all fields except `end_motifs`, `frag_len_hist`, `start_pos_hist`, `end_pos_hist`) are appended automatically. Call this directly if you want to inspect the feature matrix outside the LOO loop.

---

## Plots — `cfclassify.plots`

All plot functions return a `matplotlib.figure.Figure`. Pass `out_path` to save to disk, or omit it to get the figure object for further manipulation.

| Function | Output | Key design notes |
|---|---|---|
| `plot_end_motifs(samples, top_n, out_path)` | Grouped bar chart | Top-N motifs ranked by mean frequency across all samples; SEM error bars |
| `plot_frag_lengths(samples, out_path)` | Line plot | Normalised per-sample before averaging; reference lines at 147 bp and 200 bp |
| `plot_methylation(samples, out_path)` | Bar + scatter | Bars show group means; scatter overlays individual sample points with jitter |
| `plot_position_dist(samples, feature_key, title, out_path)` | Line plot | `feature_key` is `"start_pos_hist"` or `"end_pos_hist"`; x-axis in Mbp |

---

## Model I/O — `cfclassify.model`

### `save_model(bundle, path)` / `load_model(path)`

Serialise/deserialise a `ModelBundle` using `joblib`. The bundle is self-describing: it carries the k-mer length, contig, and feature key ordering needed to classify new samples without any external configuration.

### `save_feature_cache(samples, model_path)` / `load_feature_cache(model_path)`

Write/read the feature cache JSON at `<model_path>.features.json`. The cache stores extracted feature dicts alongside labels, so the `update` subcommand can retrain on all historical samples without re-processing BAMs.

---

## CLI entry point — `cfclassify.__main__`

Three subcommands: `train`, `predict`, `update`. See [Usage → Running](../usage.md#running) for full flag reference.
