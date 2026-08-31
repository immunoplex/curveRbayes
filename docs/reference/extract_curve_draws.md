# Posterior draw vectors for one plate's curve parameters (a,b,c,d,g)

Posterior draw vectors for one plate's curve parameters (a,b,c,d,g)

## Usage

``` r
extract_curve_draws(bayes_fit, curve_idx)
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

- curve_idx:

  Integer plate index (1-based Stan index).

## Value

Named list `term -> numeric()` (c_par emitted as c). For calib_draws
(param_scope = "curve"); emit only when persist_draws.
