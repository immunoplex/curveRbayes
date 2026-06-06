# Extract a Per-Curve Summary Table from Bayesian Results

Works with both `calibration_result_multiplate` (new format) and single
`calibration_result` (legacy format).

## Usage

``` r
summary_table_bayes(cr)
```

## Arguments

- cr:

  A `calibration_result_multiplate` or `calibration_result` with
  `method = "bayesian"`.

## Value

Data frame with one row per curve_id.
