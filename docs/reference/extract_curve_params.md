# Extract Curve-Level Posterior Summaries

Computes posterior mean, SD, and quantiles for one curve's parameters.

## Usage

``` r
extract_curve_params(bayes_fit, curve_idx = 1L, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

- curve_idx:

  Integer. Which curve (1-based Stan index).

- probs:

  Numeric vector of quantiles. Default c(0.025, 0.5, 0.975).

## Value

Data frame with columns: term, mean, sd, q2.5, q50, q97.5.
