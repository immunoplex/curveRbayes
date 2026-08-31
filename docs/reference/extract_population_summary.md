# Summarise population / noise parameters from a Bayesian fit

Summarise population / noise parameters from a Bayesian fit

## Usage

``` r
extract_population_summary(bayes_fit, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

- probs:

  Quantiles. Default c(0.025, 0.5, 0.975).

## Value

data.frame(term, estimate, std_error, q_lo, q_med, q_hi) — the
curveRcore `population$params` shape. One row per group-level scalar.
