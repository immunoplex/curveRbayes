# Per-fit sampler diagnostics from a Bayesian fit

Per-fit sampler diagnostics from a Bayesian fit

## Usage

``` r
extract_fit_diag(bayes_fit)
```

## Arguments

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

## Value

Named list matching the calib_fit_diag Bayesian columns. `converged` is
left NULL here (the caller supplies it from the ensemble's converged
flag); rhat/ess come from the CmdStanMCMC \$summary(),
divergences/treedepth/ ebfmi from bf\$diagnostics, timing from \$time(),
seed from \$metadata().
