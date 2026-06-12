# Collect All Blank Data from Bayesian Results

Extracts the per-curve blank data frames stored in each
`calibration_result$blanks` slot and stacks them into a single data
frame. Useful for QA checks on blank signal levels and for verifying
that blanks were available for every plate.

## Usage

``` r
collect_blanks_bayes(cr)
```

## Arguments

- cr:

  A `calibration_result_multiplate` or `calibration_result`.

## Value

Data frame with a `curve_id` column prepended, or NULL if no blanks are
stored in any plate.
