# Compute LOO-CV for a Fitted Bayesian Model

Extracts the log_lik generated quantity and computes PSIS-LOO.

## Usage

``` r
compute_loo(bayes_fit)
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

## Value

A `loo` object (from the loo package).
