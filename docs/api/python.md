# Python API — cfanalysis

All public symbols live in the `cfanalysis` package. Import paths are shown on each entry.

---

## Types — `cfanalysis.types`

Two `TypedDict` classes define the data schema flowing through the Python layer.

```python title="cfanalysis/types.py"
--8<-- "cfanalysis/types.py"
```

**`Features`** is `total=False`: every key is optional. Missing data (e.g. no CpG coverage in a region) is represented by the absence of the corresponding key — not by a `NaN` value. This keeps downstream code honest: a consumer must explicitly handle the absent-key case rather than silently propagating NaN.

**`Sample`** wraps a single processed BAM file.

---

## Feature conversion — `cfanalysis.features`

```python title="cfanalysis/features.py"
--8<-- "cfanalysis/features.py"
```

`metrics_to_features()` is the **sole extension point** for adding new features to the pipeline. When a new field is added to `RegionMetrics` on the C++ side:

1. Add the conversion logic here.
2. If the new field is a flat numeric scalar, `build_feature_matrix()` in `classify.py` picks it up automatically — no other changes needed.
3. If the new field is a histogram or dict (like `end_motifs`), add explicit handling to both `metrics_to_features()` and `build_feature_matrix()`.

---

## Classification — `cfanalysis.classify`

```python title="cfanalysis/classify.py"
--8<-- "cfanalysis/classify.py"
```

### `run_loo_cv(samples, top_n=20)`

Runs leave-one-out cross-validation and returns a dict:

```python
{
    "predictions": list[str],   # one label per sample, in input order
    "labels":      list[str],   # ground-truth labels in the same order
    "report":      str,         # sklearn classification_report string
}
```

The three leakage guards (motif selection, scaling, hyperparameter search) are all applied inside the fold loop. See [Classification — Implementation](../classification.md#implementation) for an annotated walkthrough.

### `build_feature_matrix(samples, motifs)`

Builds `(X, y, feature_names)` from a sample list and an explicit motif list. Flat numeric keys in `Features` (all fields except `end_motifs`, `frag_len_hist`, `start_pos_hist`, `end_pos_hist`) are appended automatically. Call this directly if you want to inspect the feature matrix outside the LOO loop.

---

## Plots — `cfanalysis.plots`

All plot functions return a `matplotlib.figure.Figure`. Pass `out_path` to save to disk, or omit it to get the figure object for further manipulation.

```python title="cfanalysis/plots.py"
--8<-- "cfanalysis/plots.py"
```

| Function | Output | Key design notes |
|---|---|---|
| `plot_end_motifs(samples, top_n, out_path)` | Grouped bar chart | Top-N motifs ranked by mean frequency across all samples; SEM error bars |
| `plot_frag_lengths(samples, out_path)` | Line plot | Normalised per-sample before averaging; reference lines at 147 bp and 200 bp |
| `plot_methylation(samples, out_path)` | Bar + scatter | Bars show group means; scatter overlays individual sample points with jitter |
| `plot_position_dist(samples, feature_key, title, out_path)` | Line plot | `feature_key` is `"start_pos_hist"` or `"end_pos_hist"`; x-axis in Mbp |

---

## CLI entry point — `cfanalysis.__main__`

```python title="cfanalysis/__main__.py"
--8<-- "cfanalysis/__main__.py"
```

The manifest CSV is read with `csv.DictReader` — any extra columns are ignored. BAM paths relative to the manifest directory are resolved before being passed to `cfextract.extract_features()`.
