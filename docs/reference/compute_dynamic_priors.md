# Compute Data-Adaptive Priors for Stan Models

Inspects preprocessed standard curve data and produces all hyperprior
locations and scales needed by the Stan model `data {}` block.

## Usage

``` r
compute_dynamic_priors(
  data,
  response_variable,
  fixed_a = NULL,
  model_family = "logistic4"
)
```

## Arguments

- data:

  Data frame of preprocessed standards with `concentration` and response
  columns (already on the fitting scale).

- response_variable:

  Character. Response column name.

- fixed_a:

  Numeric or NULL. If non-NULL, the prior on `a` is tightened to a soft
  constraint centered on this value (on the fitting scale).

- model_family:

  Character. One of the curveRcore model names. Default `"logistic4"`.

## Value

Named list of hyperprior values matching the Stan `data` block.
