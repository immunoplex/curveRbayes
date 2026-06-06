# Compare Models via LOO-CV and Stacking Weights

Given a list of fitted Bayesian models (one per model family), computes
LOO for each and Bayesian stacking weights.

## Usage

``` r
compare_models_loo(fits)
```

## Arguments

- fits:

  Named list of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md)
  outputs.

## Value

Named list with `best_model_name`, `criterion`, `loo_results`,
`comparison`, `weights`.
