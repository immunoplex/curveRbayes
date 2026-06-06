# Fit a Single Model Family via MCMC

Runs HMC/NUTS sampling for one model family across all curve_ids in the
data. Returns a structured list containing the CmdStanMCMC fit,
posterior draws, and metadata.

## Usage

``` r
fit_bayes_single(
  stan_data,
  model_family = "logistic4",
  chains = 4L,
  warmup = 1000L,
  sampling = 1000L,
  adapt_delta = 0.9,
  max_treedepth = 12L,
  seed = NULL,
  compiled_model = NULL,
  verbose = FALSE
)
```

## Arguments

- stan_data:

  Named list from
  [`build_stan_data()`](https://immunoplex.github.io/curveRbayes/reference/build_stan_data.md).

- model_family:

  Character. Model name.

- chains:

  Integer. Default 4.

- warmup:

  Integer. Default 1000.

- sampling:

  Integer. Default 1000.

- adapt_delta:

  Numeric. Default 0.9.

- max_treedepth:

  Integer. Default 12.

- seed:

  Integer or NULL.

- compiled_model:

  Optional pre-compiled CmdStanModel.

- verbose:

  Logical.

## Value

A named list with `model_family`, `fit`, `draws`, `n_curves`,
`stan_data`, `diagnostics`.
