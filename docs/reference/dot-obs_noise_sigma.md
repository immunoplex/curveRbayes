# Observation-noise sigma for one posterior draw (shared by grid + samples)

Single definition of the per-draw observation SD so the grid and the
sample back-calculation can never diverge. Homoscedastic returns the
shared `sigma_obs` draw; heteroscedastic returns the O'Malley
power-of-mean form `exp(log_sigma0 + log_sigma_slope * log(|mu|))`.

## Usage

``` r
.obs_noise_sigma(
  mu,
  s,
  use_hetero,
  sigma_obs_draws,
  log_sigma0_draws = NULL,
  log_sigma_slope_draws = NULL
)
```

## Arguments

- mu:

  Numeric. Mean response the noise is centred on (grid: the forward mean
  at the grid point; samples: the observed response, the best available
  proxy for the underlying mean).

- s:

  Integer. Draw index.

- use_hetero:

  Logical. Whether the heteroscedastic model was fitted.

- sigma_obs_draws, log_sigma0_draws, log_sigma_slope_draws:

  Numeric draw vectors (the last two may be NULL under the homoscedastic
  model).
