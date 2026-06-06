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

- grainsize:

  Integer. reduce_sum grain size. Default 1.

## Value

Named list suitable for `cmdstanr::sample(data = ...)`.
