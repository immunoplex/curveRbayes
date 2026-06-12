# Predict Grid Response from Posterior Draws (Bayesian)

For each grid point, evaluates the forward model at every posterior
draw, adds observation noise, then back-calculates concentration to
produce a precision profile.

## Usage

``` r
predict_grid_bayes(
  grid,
  bayes_fit,
  curve_idx = 1L,
  n_draws = NULL,
  cv_x_max = 150,
  pcov_threshold = 20,
  is_log_x = TRUE,
  is_log_response = TRUE
)
```

## Arguments

- grid:

  Data frame from
  [`curveRcore::generate_prediction_grid()`](https://immunoplex.github.io/curveRcore/reference/generate_prediction_grid.html).

- bayes_fit:

  Output of
  [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md).

- curve_idx:

  Integer. Which curve (1-based Stan index).

- n_draws:

  Integer or NULL. Subsample this many draws.

- cv_x_max:

  Numeric. Cap for pcov/pcov_rmse. Default 150.

- pcov_threshold:

  Numeric. Percent CV threshold for pcov_pass. Default 20.

- is_log_x:

  Logical. Default TRUE.

- is_log_response:

  Logical. Whether the response is log10-transformed. Passed to
  [`curveRcore::enrich_grid_with_d2y()`](https://immunoplex.github.io/curveRcore/reference/enrich_grid_with_d2y.html)
  for second-derivative enrichment. Default TRUE.

## Value

`grid` with added columns: `predicted_response`, `ci_lower`, `ci_upper`,
`predicted_concentration`, `se_concentration`, `pcov`, `pcov_rmse`,
`pcov_pass`, `noise_mode`.

## Details

When the model was fitted with `use_heteroscedastic_noise = TRUE`, the
noise injected at Step 2 scales with the predicted response magnitude
(`sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))`), giving
the O'Malley (2008) CDAN precision profile. When
`use_heteroscedastic_noise = FALSE`, a constant `sigma_obs` is used and
the profile reflects posterior-predictive uncertainty.
