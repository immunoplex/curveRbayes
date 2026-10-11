# Back-Calculate Sample Concentrations from Posterior Draws

For each test sample, evaluates the inverse model at every posterior
draw to produce a posterior distribution of predicted concentration. The
variance definition matches
[`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
via `include_measurement_error`, so the sample pcov is directly
comparable to the precision profile.

## Usage

``` r
predict_samples_bayes(
  samples,
  bayes_fit,
  curve_idx = 1L,
  response_variable,
  is_log_response = TRUE,
  n_draws = NULL,
  cv_x_max = 150,
  pcov_threshold = 20,
  is_log_x = TRUE,
  include_measurement_error = TRUE,
  lloq_log10 = NA_real_,
  uloq_log10 = NA_real_,
  decision_cutoff = NULL
)
```

## Arguments

- samples:

  Data frame of test samples.

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

- curve_idx:

  Integer. Which curve (1-based Stan index).

- response_variable:

  Character.

- is_log_response:

  Logical.

- n_draws:

  Integer or NULL.

- cv_x_max:

  Numeric. Default 150.

- pcov_threshold:

  Numeric. Percent CV threshold for pcov_pass. Default 20.

- is_log_x:

  Logical. Default TRUE.

- include_measurement_error:

  Logical. If TRUE (default) inject the SAME observation noise the grid
  uses before back-calculating, so samples lie on the measurement/CDAN
  profile. If FALSE, invert the observed response with no added noise
  (curve/parameter precision only). MUST match the value passed to
  [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
  for the same fit.

- lloq_log10, uloq_log10:

  Numeric. This curve's LLOQ/ULOQ (log10, pre-dilution) for
  `p_below_lloq` / `p_above_uloq`. Default `NA`.

- decision_cutoff:

  Numeric or NULL. Dilution-corrected concentration for
  `p_below_cutoff`. Default NULL (column is `NA`).

## Value

Data frame with original sample columns plus prediction columns,
including the predictive interval / probability columns from
[`curveRcore::summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.html).
Those are computed from ALL draws, with draws whose response falls
beyond an asymptote kept as censored; `se_concentration` / `pcov` are
unchanged and still use the invertible draws only.
