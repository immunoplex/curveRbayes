# =============================================================================
# predict_bayes.R — Posterior predictive grid and sample back-calculation
#
# Implements two precision profile modes selectable via the
# `use_heteroscedastic_noise` argument of fit_calibration_bayes():
#
# Mode 0 — Posterior-predictive precision (homoscedastic):
#   sigma_i is constant (sigma_obs from the posterior). The precision
#   profile reflects mainly the inverse-curve geometry.
#
# Mode 1 — CDAN precision (heteroscedastic, O'Malley 2008):
#   sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))
#   A power-of-mean variance function where noise scales with signal level.
#   For each posterior draw, the grid point, and each CDAN replicate:
#     1. Evaluate the forward model  →  mu_s
#     2. Compute sigma_i from the heteroscedastic noise model
#     3. Draw  y_noisy ~ Student-t(nu, mu_s, sigma_i)
#     4. Back-calculate concentration from y_noisy  →  x_hat
#   pcov_rmse is then the relative RMSE of x_hat around x_true, which is
#   the O'Malley CDAN precision measure.
# =============================================================================


#' Predict Grid Response from Posterior Draws (Bayesian)
#'
#' For each grid point, evaluates the forward model at every posterior
#' draw, adds observation noise, then back-calculates concentration to
#' produce a precision profile.
#'
#' When the model was fitted with `use_heteroscedastic_noise = TRUE`,
#' the noise injected at Step 2 scales with the predicted response
#' magnitude (`sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))`),
#' giving the O'Malley (2008) CDAN precision profile. When
#' `use_heteroscedastic_noise = FALSE`, a constant `sigma_obs` is used
#' and the profile reflects posterior-predictive uncertainty.
#'
#' @param grid Data frame from [curveRcore::generate_prediction_grid()].
#' @param bayes_fit Output of [curveRbayes::fit_bayes_single()].
#' @param curve_idx Integer. Which curve (1-based Stan index).
#' @param n_draws Integer or NULL. Subsample this many draws.
#' @param cv_x_max Numeric. Cap for pcov/pcov_rmse. Default 150.
#' @param pcov_threshold Numeric. Percent CV threshold for pcov_pass. Default 20.
#' @param is_log_x Logical. Default TRUE.
#' @param is_log_response Logical. Whether the response is log10-transformed.
#'   Passed to [curveRcore::enrich_grid_with_d2y()] for second-derivative
#'   enrichment. Default TRUE.
#'
#' @return `grid` with added columns: `predicted_response`, `ci_lower`,
#'   `ci_upper`, `predicted_concentration`, `se_concentration`, `pcov`,
#'   `pcov_rmse`, `pcov_pass`, `noise_mode`.
#'
#' @export
predict_grid_bayes <- function(grid, bayes_fit, curve_idx = 1L,
                               n_draws = NULL, cv_x_max = 150,
                               pcov_threshold = 20,
                               is_log_x = TRUE,
                               is_log_response = TRUE
                               ) {

  family <- bayes_fit$model_family
  draws  <- bayes_fit$draws
  p      <- curve_idx
  n_grid <- nrow(grid)

  # ── Curve shape draws ──
  a_draws <- as.numeric(draws[[paste0("a[", p, "]")]])
  b_draws <- as.numeric(draws[[paste0("b[", p, "]")]])
  c_draws <- as.numeric(draws[[paste0("c_par[", p, "]")]])
  d_draws <- as.numeric(draws[[paste0("d[", p, "]")]])
  g_draws <- if (family %in% c("logistic5", "loglogistic5"))
    as.numeric(draws[[paste0("g[", p, "]")]]) else NULL

  # ── Noise draws — homoscedastic path ──
  sigma_obs_draws <- as.numeric(draws[["sigma_obs"]])
  nu_draws        <- as.numeric(draws[["nu"]])

  # ── Noise draws — heteroscedastic (CDAN) path ──
  # These columns are present whenever use_heteroscedastic_noise was set
  # in build_stan_data() (the Stan model always estimates them).
  has_hetero      <- all(c("log_sigma0", "log_sigma_slope") %in% names(draws))
  log_sigma0_draws     <- if (has_hetero) as.numeric(draws[["log_sigma0"]])     else NULL
  log_sigma_slope_draws <- if (has_hetero) as.numeric(draws[["log_sigma_slope"]]) else NULL

  # Determine which noise path was active at fitting time.
  # bayes_fit$stan_data$use_heteroscedastic_noise is 0 or 1.
  use_hetero <- isTRUE(bayes_fit$stan_data$use_heteroscedastic_noise == 1L)

  S <- length(a_draws)
  if (!is.null(n_draws) && n_draws < S) {
    idx <- sample.int(S, n_draws)
    a_draws <- a_draws[idx]; b_draws <- b_draws[idx]
    c_draws <- c_draws[idx]; d_draws <- d_draws[idx]
    if (!is.null(g_draws))            g_draws            <- g_draws[idx]
    sigma_obs_draws <- sigma_obs_draws[idx]
    nu_draws        <- nu_draws[idx]
    if (!is.null(log_sigma0_draws)) {
      log_sigma0_draws      <- log_sigma0_draws[idx]
      log_sigma_slope_draws <- log_sigma_slope_draws[idx]
    }
    S <- n_draws
  }

  # Forward model dispatcher
  fwd <- switch(family,
                logistic4    = function(x, s) curveRcore::logistic4(x, a_draws[s], b_draws[s],
                                                                    c_draws[s], d_draws[s]),
                logistic5    = function(x, s) curveRcore::logistic5(x, a_draws[s], b_draws[s],
                                                                    c_draws[s], d_draws[s], g_draws[s]),
                loglogistic4 = function(x, s) curveRcore::loglogistic4(x, a_draws[s], b_draws[s],
                                                                       c_draws[s], d_draws[s]),
                loglogistic5 = function(x, s) curveRcore::loglogistic5(x, a_draws[s], b_draws[s],
                                                                       c_draws[s], d_draws[s], g_draws[s]),
                gompertz4    = function(x, s) curveRcore::gompertz4(x, a_draws[s], b_draws[s],
                                                                    c_draws[s], d_draws[s])
  )

  # Inverse model dispatcher
  inv <- switch(family,
                logistic4    = function(y, s) curveRcore::inv_logistic4(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s]),
                logistic5    = function(y, s) curveRcore::inv_logistic5(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s], g_draws[s]),
                loglogistic4 = function(y, s) curveRcore::inv_loglogistic4(y, a_draws[s], b_draws[s],
                                                                           c_draws[s], d_draws[s]),
                loglogistic5 = function(y, s) curveRcore::inv_loglogistic5(y, a_draws[s], b_draws[s],
                                                                           c_draws[s], d_draws[s], g_draws[s]),
                gompertz4    = function(y, s) curveRcore::inv_gompertz4(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s])
  )

  # Pre-allocate matrices
  y_mat <- matrix(NA_real_, nrow = n_grid, ncol = S)
  x_mat <- matrix(NA_real_, nrow = n_grid, ncol = S)

  for (s in seq_len(S)) {
    for (i in seq_len(n_grid)) {
      # Step 1: Forward prediction (exact mean)
      mu_s <- suppressWarnings(
        tryCatch(fwd(grid$x_fit[i], s), error = function(e) NA_real_)
      )
      y_mat[i, s] <- mu_s

      if (!is.finite(mu_s)) {
        x_mat[i, s] <- NA_real_
        next
      }

      # Step 2: Draw noisy observation.
      # Heteroscedastic (CDAN): sigma scales with |mu_s|.
      # Homoscedastic:          sigma is the shared constant sigma_obs.
      sigma_s <- if (use_hetero && !is.null(log_sigma0_draws)) {
        exp(log_sigma0_draws[s] + log_sigma_slope_draws[s] * log(abs(mu_s) + 1e-10))
      } else {
        sigma_obs_draws[s]
      }
      y_noisy <- mu_s + sigma_s * stats::rt(1, df = nu_draws[s])

      # Step 3: Back-calculate from noisy response
      x_mat[i, s] <- suppressWarnings(
        tryCatch(inv(y_noisy, s), error = function(e) NA_real_)
      )
    }
  }

  # ── Summaries ──
  grid$predicted_response      <- rowMeans(y_mat, na.rm = TRUE)
  grid$ci_lower                <- apply(y_mat, 1, stats::quantile,
                                        probs = 0.025, na.rm = TRUE)
  grid$ci_upper                <- apply(y_mat, 1, stats::quantile,
                                        probs = 0.975, na.rm = TRUE)
  grid$predicted_concentration <- apply(x_mat, 1, stats::median, na.rm = TRUE)
  grid$se_concentration        <- apply(x_mat, 1, stats::sd, na.rm = TRUE)

  # pcov: posterior SD of back-calculated concentration as CV (%)
  grid$pcov <- vapply(seq_len(n_grid), function(i) {
    se_i <- grid$se_concentration[i]
    if (!is.finite(se_i)) return(cv_x_max)
    raw_cv <- if (is_log_x) se_i * log(10) * 100
    else if (abs(grid$predicted_concentration[i]) > 1e-10)
      (se_i / abs(grid$predicted_concentration[i])) * 100
    else Inf
    min(raw_cv, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  # pcov_rmse: relative RMSE of back-calculated concentration around x_true.
  # When use_heteroscedastic_noise = TRUE, this is the O'Malley (2008) CDAN
  # precision measure — noise scales with signal level, so the profile
  # captures both inverse-curve geometry and concentration-dependent noise.
  # When use_heteroscedastic_noise = FALSE, this is still a valid
  # posterior-predictive precision summary but does not constitute CDAN.
  grid$pcov_rmse <- vapply(seq_len(n_grid), function(i) {
    x_draws <- x_mat[i, ]
    x_true  <- grid$x_fit[i]
    ok <- is.finite(x_draws)
    if (sum(ok) < 2) return(cv_x_max)
    mse <- mean((x_draws[ok] - x_true)^2)
    rmse <- sqrt(mse)
    raw_rrmse <- if (is_log_x) rmse * log(10) * 100
    else if (abs(x_true) > 1e-10) (rmse / abs(x_true)) * 100
    else Inf
    min(raw_rrmse, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  # pcov_pass uses pcov_threshold (precision budget), NOT cv_x_max (hard cap).
  grid$pcov_pass <- !is.na(grid$pcov) & grid$pcov < pcov_threshold

  # Record which noise mode was used (useful for plotting / auditing)
  grid$noise_mode <- if (use_hetero) "heteroscedastic" else "homoscedastic"

  grid <- curveRcore::enrich_grid_with_d2y(grid, is_log_response = is_log_response)

  grid
}


#' Back-Calculate Sample Concentrations from Posterior Draws
#'
#' For each test sample, evaluates the inverse model at every posterior
#' draw to produce a full posterior distribution of predicted concentration.
#'
#' @param samples Data frame of test samples.
#' @param bayes_fit Output of [curveRbayes::fit_bayes_single()].
#' @param curve_idx Integer. Which curve (1-based Stan index).
#' @param response_variable Character.
#' @param is_log_response Logical.
#' @param n_draws Integer or NULL.
#' @param cv_x_max Numeric. Default 150.
#' @param pcov_threshold Numeric. Percent CV threshold for pcov_pass. Default 20.
#' @param is_log_x Logical. Default TRUE.
#'
#' @return Data frame with original sample columns plus prediction columns.
#'
#' @export
predict_samples_bayes <- function(samples, bayes_fit, curve_idx = 1L,
                                  response_variable,
                                  is_log_response = TRUE,
                                  n_draws = NULL, cv_x_max = 150,
                                  pcov_threshold = 20,
                                  is_log_x = TRUE) {

  family <- bayes_fit$model_family
  draws  <- bayes_fit$draws
  p      <- curve_idx
  n      <- nrow(samples)
  if (n == 0) return(samples)

  a_draws <- as.numeric(draws[[paste0("a[", p, "]")]])
  b_draws <- as.numeric(draws[[paste0("b[", p, "]")]])
  c_draws <- as.numeric(draws[[paste0("c_par[", p, "]")]])
  d_draws <- as.numeric(draws[[paste0("d[", p, "]")]])
  g_draws <- if (family %in% c("logistic5", "loglogistic5"))
    as.numeric(draws[[paste0("g[", p, "]")]]) else NULL

  S <- length(a_draws)
  if (!is.null(n_draws) && n_draws < S) {
    idx <- sample.int(S, n_draws)
    a_draws <- a_draws[idx]; b_draws <- b_draws[idx]
    c_draws <- c_draws[idx]; d_draws <- d_draws[idx]
    if (!is.null(g_draws)) g_draws <- g_draws[idx]
    S <- n_draws
  }

  inv <- switch(family,
                logistic4    = function(y, s) curveRcore::inv_logistic4(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s]),
                logistic5    = function(y, s) curveRcore::inv_logistic5(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s], g_draws[s]),
                loglogistic4 = function(y, s) curveRcore::inv_loglogistic4(y, a_draws[s], b_draws[s],
                                                                           c_draws[s], d_draws[s]),
                loglogistic5 = function(y, s) curveRcore::inv_loglogistic5(y, a_draws[s], b_draws[s],
                                                                           c_draws[s], d_draws[s], g_draws[s]),
                gompertz4    = function(y, s) curveRcore::inv_gompertz4(y, a_draws[s], b_draws[s],
                                                                        c_draws[s], d_draws[s])
  )

  raw_response <- samples[[response_variable]]
  fit_response <- if (is_log_response) log10(pmax(raw_response, 1e-6)) else raw_response

  # Posterior matrix: rows = samples, cols = draws
  # For samples, the response IS the noisy observation — no need to add noise
  x_mat <- matrix(NA_real_, nrow = n, ncol = S)
  for (s in seq_len(S)) {
    for (i in seq_len(n)) {
      x_mat[i, s] <- suppressWarnings(
        tryCatch(inv(fit_response[i], s), error = function(e) NA_real_)
      )
    }
  }

  x_median <- apply(x_mat, 1, stats::median, na.rm = TRUE)
  x_se     <- apply(x_mat, 1, stats::sd, na.rm = TRUE)

  dilution <- if ("dilution" %in% names(samples)) samples$dilution else 1

  samples$raw_assay_response            <- raw_response
  samples$observed_response_fit         <- fit_response
  samples$predicted_log10_concentration <- if (is_log_x) x_median else log10(pmax(x_median, 1e-20))
  samples$predicted_concentration       <- x_median
  samples$final_concentration           <- if (is_log_x) 10^x_median * dilution else x_median * dilution
  samples$se_concentration              <- x_se

  # pcov: CV (%)
  samples$pcov <- vapply(seq_len(n), function(i) {
    if (!is.finite(x_se[i])) return(cv_x_max)
    raw_cv <- if (is_log_x) x_se[i] * log(10) * 100
    else if (abs(x_median[i]) > 1e-10) (x_se[i] / abs(x_median[i])) * 100
    else Inf
    min(raw_cv, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  # pcov_rmse: for samples, no known truth, so RMSE = SE (bias undefined)
  samples$pcov_rmse <- samples$pcov

  # pcov_pass uses pcov_threshold (precision budget), NOT cv_x_max (hard cap).
  # cv_x_max is an overflow guard; marking everything below 150 % as "pass"
  # would be misleading — the usable range is defined by pcov_threshold.
  samples$pcov_pass <- !is.na(samples$pcov) & samples$pcov < pcov_threshold
  samples
}
