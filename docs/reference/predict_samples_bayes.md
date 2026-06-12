# Back-Calculate Sample Concentrations from Posterior Draws

For each test sample, evaluates the inverse model at every posterior
draw to produce a full posterior distribution of predicted
concentration.

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
  is_log_x = TRUE
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

## Value

Data frame with original sample columns plus prediction columns.
