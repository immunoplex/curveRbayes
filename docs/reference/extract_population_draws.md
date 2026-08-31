# Posterior draw vectors for the population / noise scalars

Posterior draw vectors for the population / noise scalars

## Usage

``` r
extract_population_draws(bayes_fit)
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

## Value

Named list `term -> numeric()`, iteration-ordered (all terms share one
order, so they can be column-bound into a joint posterior). For
calib_draws (param_scope = "population"); emit only when persist_draws.
