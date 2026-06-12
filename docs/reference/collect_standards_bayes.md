# Collect All Standard Data from Bayesian Results

Extracts the per-curve standards data frames stored in each
`calibration_result$standards` slot and stacks them into a single data
frame. Useful for verifying standard coverage or plotting observed data
alongside the fitted curve.

## Usage

``` r
collect_standards_bayes(cr)
```

## Arguments

- cr:

  A `calibration_result_multiplate` or `calibration_result`.

## Value

Data frame with a `curve_id` column prepended, or NULL if no standards
are stored in any plate.
