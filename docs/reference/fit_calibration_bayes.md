# Fit Bayesian Hierarchical Calibration Curves

Fits a hierarchical Bayesian model to preprocessed standard curve data
across one or more curve_ids simultaneously. Returns a
`calibration_result_multiplate` with one entry per curve_id, each
containing its own grid predictions, sample predictions, and parameter
summaries.

## Usage

``` r
fit_calibration_bayes(
  standards,
  samples = NULL,
  blanks = NULL,
  response_var,
  model_names = c("logistic4", "gompertz4"),
  is_log_response = TRUE,
  is_log_independent = TRUE,
  std_curve_conc,
  fixed_a = NULL,
  cv_x_max = 150,
  pcov_threshold = 20,
  min_dynamic_range_log10 = 0.5,
  max_rel_se = 5,
  n_grid = 200L,
  grid_min_conc = 1e-04,
  grid_max_conc = NULL,
  chains = 4L,
  warmup = 1000L,
  sampling = 1000L,
  adapt_delta = 0.9,
  seed = NULL,
  n_draws_predict = 500L,
  n_draws_ensemble = 260L,
  compute_all_grids = FALSE,
  use_heteroscedastic_noise = FALSE,
  run_loo = NULL,
  verbose = FALSE
)
```

## Arguments

- standards:

  Data frame. Preprocessed stacked standard curve data. Must contain
  `curve_id`, a response column, and a `concentration` column — all on
  the fitting scale.

- samples:

  Data frame or NULL. Stacked sample data with `curve_id` and the
  response column (on the raw measurement scale).

- blanks:

  Data frame or NULL. Blank well data with `curve_id` and the response
  column (on the fitting scale). When supplied, blanks are passed to
  Stan to anchor the lower asymptote via a separate likelihood term, and
  are stored in each per-curve `calibration_result$blanks` slot for
  downstream QA. Default NULL.

- response_var:

  Character. Name of the response column.

- model_names:

  Character vector. Models to fit. Default
  `c("logistic4", "gompertz4")`.

- is_log_response:

  Logical. Default TRUE.

- is_log_independent:

  Logical. Default TRUE.

- std_curve_conc:

  Numeric. Undiluted standard concentration.

- fixed_a:

  Numeric or NULL. Fixed lower asymptote (fitting scale).

- cv_x_max:

  Numeric. Default 150.

- pcov_threshold:

  Numeric. Percent CV threshold for LLOQ/ULOQ determination and the
  dynamic-range eligibility gate. Default 20.

- min_dynamic_range_log10:

  Numeric. Minimum dynamic range (log10) for eligibility. Default 0.5.

- max_rel_se:

  Numeric. Maximum relative SE (SD/\|mean\|) permitted for any
  parameter. Default 5.0.

- n_grid:

  Integer. Default 200.

- grid_min_conc:

  Numeric. Default 1e-4.

- grid_max_conc:

  Numeric or NULL.

- chains:

  Integer. Default 4.

- warmup:

  Integer. Default 1000.

- sampling:

  Integer. Default 1000.

- adapt_delta:

  Numeric. Default 0.9.

- seed:

  Integer or NULL.

- n_draws_predict:

  Integer. Number of posterior draws for the best-model grid and sample
  predictions. Default 500.

- n_draws_ensemble:

  Integer. Number of posterior draws for non-best-model precision grids.
  Default 260.

- compute_all_grids:

  Logical. If TRUE, compute precision grids for every converged model.
  Required for eligibility gating when more than one model is fitted.
  Default FALSE.

- use_heteroscedastic_noise:

  Logical. If TRUE, the Stan models use a power-of-mean variance
  function `sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))`
  in the likelihood, and the same sigma_i is injected when generating
  the CDAN noisy observations in
  [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md).
  This restores the O'Malley (2008) CDAN precision profile
  interpretation. If FALSE (default), a constant `sigma_obs` is used and
  the precision profiles reflect posterior-predictive uncertainty driven
  mainly by inverse- curve geometry.

- run_loo:

  Logical or NULL. Default NULL (auto).

- verbose:

  Logical. Default FALSE.

## Value

A `calibration_result_multiplate` object (from curveRcore). Each
per-plate `$selection` contains `$assessments`, `$eligible_models`, and
`$fallback` from the eligibility gating. Each per-plate
`calibration_result` carries `$standards` and `$blanks` slots with the
per-curve subsets of the input data.
