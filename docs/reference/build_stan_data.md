# Build Stan Data List for a Model Family

Converts a preprocessed stacked data frame into the named list format
expected by the Stan model's `data {}` block.

## Usage

``` r
build_stan_data(
  standards,
  response_variable,
  priors,
  model_family = "logistic4",
  curve_id_map,
  blanks = NULL,
  use_heteroscedastic_noise = FALSE,
  grainsize = 1L
)
```

## Arguments

- standards:

  Data frame. Preprocessed stacked standards with `curve_id`,
  `concentration`, and response column.

- response_variable:

  Character. Response column name.

- priors:

  Named list from
  [`compute_dynamic_priors()`](https://immunoplex.github.io/curveRbayes/reference/compute_dynamic_priors.md).

- model_family:

  Character. One of the curveRcore model names.

- curve_id_map:

  Named integer vector mapping curve_id values to 1-based Stan indices.

- blanks:

  Data frame or NULL. Blank well data. Must contain `curve_id` and the
  response column (on the fitting scale). Passed directly to Stan to
  anchor the lower asymptote.

- use_heteroscedastic_noise:

  Logical. If TRUE, passes `use_heteroscedastic_noise = 1L` to Stan,
  enabling the power-of-mean noise model
  `sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))`. If FALSE
  (default), the homoscedastic `sigma_obs` path is used.

- grainsize:

  Integer. reduce_sum grain size. Default 1. NOTE: grainsize is declared
  in the Stan data block but reduce_sum is not yet active in the model
  bodies; this field is reserved for a future threading implementation
  and is currently ignored by Stan.

## Value

Named list suitable for `cmdstanr::sample(data = ...)`.
