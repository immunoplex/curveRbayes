# =============================================================================
# predict_bayes.R — Posterior predictive grid and sample back-calculation
#
# ── PRECISION DEFINITION (single source of truth) ────────────────────────────
# The precision profile (grid) and the per-sample pcov MUST be computed under
# the SAME error definition, or sample points will not lie on the profile. That
# consistency is now controlled by ONE switch, `include_measurement_error`,
# which is threaded identically into BOTH predict_grid_bayes() and
# predict_samples_bayes():
#
#   include_measurement_error = TRUE  (DEFAULT) — "measurement" precision:
#     A fresh observation-noise draw is injected before back-calculation, so the
#     profile answers "given the assay's measurement noise, how precisely can a
#     concentration be recovered from a single reading?" This is the O'Malley
#     (2008) CDAN-style precision profile the app has always plotted. Samples get
#     the SAME noise injection, so the sample cloud lies on the profile.
#
#   include_measurement_error = FALSE — "curve" precision:
#     No observation noise is added. Precision reflects calibration-curve
#     (posterior parameter) uncertainty ONLY. Both paths invert a FIXED reference
#     response across the posterior, so the grid is non-degenerate and the sample
#     cloud again lies on the profile — just lower. Use this when the plate has
#     too few standards/controls to estimate measurement error reliably (see
#     UNDERSTANDING_precision_and_measurement_error.md).
#
# The noise MODEL used when include_measurement_error = TRUE is still selected by
# use_heteroscedastic_noise (homoscedastic sigma_obs vs the power-of-mean CDAN
# form). That switch is orthogonal and only matters when noise is included.
#
# Because both paths call the SAME .obs_noise_sigma() helper, the two can never
# drift apart again.
# =============================================================================


#' Observation-noise sigma for one posterior draw (shared by grid + samples)
#'
#' Single definition of the per-draw observation SD so the grid and the sample
#' back-calculation can never diverge. Homoscedastic returns the shared
#' `sigma_obs` draw; heteroscedastic returns the O'Malley power-of-mean form
#' `exp(log_sigma0 + log_sigma_slope * log(|mu|))`.
#'
#' @param mu Numeric. Mean response the noise is centred on (grid: the forward
#'   mean at the grid point; samples: the observed response, the best available
#'   proxy for the underlying mean).
#' @param s Integer. Draw index.
#' @param use_hetero Logical. Whether the heteroscedastic model was fitted.
#' @param sigma_obs_draws,log_sigma0_draws,log_sigma_slope_draws Numeric draw
#'   vectors (the last two may be NULL under the homoscedastic model).
#' @keywords internal
.obs_noise_sigma <- function(mu, s, use_hetero,
                             sigma_obs_draws,
                             log_sigma0_draws = NULL,
                             log_sigma_slope_draws = NULL) {
  if (isTRUE(use_hetero) && !is.null(log_sigma0_draws)) {
    exp(log_sigma0_draws[s] + log_sigma_slope_draws[s] * log(abs(mu) + 1e-10))
  } else {
    sigma_obs_draws[s]
  }
}


#' Predict Grid Response from Posterior Draws (Bayesian)
#'
#' For each grid point, evaluates the forward model at every posterior draw and
#' back-calculates concentration to produce a precision profile. Whether a fresh
#' observation-noise draw is injected before back-calculation is governed by
#' `include_measurement_error` (see file header).
#'
#' @param grid Data frame from [curveRcore::generate_prediction_grid()].
#' @param bayes_fit Output of [curveRbayes::fit_bayes_single()].
#' @param curve_idx Integer. Which curve (1-based Stan index).
#' @param n_draws Integer or NULL. Subsample this many draws.
#' @param cv_x_max Numeric. Cap for pcov/pcov_rmse. Default 150.
#' @param pcov_threshold Numeric. Percent CV threshold for pcov_pass. Default 20.
#' @param is_log_x Logical. Default TRUE.
#' @param is_log_response Logical. Whether the response is log10-transformed.
#'   Passed to [curveRcore::enrich_grid_with_d2y()]. Default TRUE.
#' @param include_measurement_error Logical. If TRUE (default) inject observation
#'   noise before back-calculation (measurement/CDAN precision). If FALSE, invert
#'   the fixed posterior-mean reference response across draws (curve/parameter
#'   precision only). See file header.
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
                               is_log_response = TRUE,
                               include_measurement_error = TRUE) {

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

  # ── Noise draws ──
  sigma_obs_draws <- as.numeric(draws[["sigma_obs"]])
  nu_draws        <- as.numeric(draws[["nu"]])
  has_hetero      <- all(c("log_sigma0", "log_sigma_slope") %in% names(draws))
  log_sigma0_draws      <- if (has_hetero) as.numeric(draws[["log_sigma0"]])      else NULL
  log_sigma_slope_draws <- if (has_hetero) as.numeric(draws[["log_sigma_slope"]]) else NULL
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

  # ── Pass 1: clean forward means (needed for predicted_response, CI, and — in
  #    curve-only mode — the fixed reference response that is inverted). ──
  y_mat <- matrix(NA_real_, nrow = n_grid, ncol = S)
  for (s in seq_len(S)) for (i in seq_len(n_grid)) {
    y_mat[i, s] <- suppressWarnings(tryCatch(fwd(grid$x_fit[i], s),
                                             error = function(e) NA_real_))
  }
  predicted_response <- rowMeans(y_mat, na.rm = TRUE)

  # ── Pass 2: back-calculate under the selected precision definition. ──
  x_mat <- matrix(NA_real_, nrow = n_grid, ncol = S)
  if (isTRUE(include_measurement_error)) {
    # Measurement / CDAN: perturb each draw's forward mean, then invert.
    for (s in seq_len(S)) for (i in seq_len(n_grid)) {
      mu_s <- y_mat[i, s]
      if (!is.finite(mu_s)) next
      sigma_s <- .obs_noise_sigma(mu_s, s, use_hetero,
                                  sigma_obs_draws, log_sigma0_draws, log_sigma_slope_draws)
      y_noisy <- mu_s + sigma_s * stats::rt(1, df = nu_draws[s])
      x_mat[i, s] <- suppressWarnings(tryCatch(inv(y_noisy, s), error = function(e) NA_real_))
    }
  } else {
    # Curve only: invert a FIXED reference response (the posterior-mean forward
    # value at the grid point) across the posterior. This isolates parameter
    # uncertainty and — unlike inverting each draw's own mean, which would round-
    # trip back to x exactly — yields a non-degenerate spread. It mirrors how the
    # sample path inverts a fixed observed response, so samples lie on the profile.
    for (s in seq_len(S)) for (i in seq_len(n_grid)) {
      y_ref <- predicted_response[i]
      if (!is.finite(y_ref)) next
      x_mat[i, s] <- suppressWarnings(tryCatch(inv(y_ref, s), error = function(e) NA_real_))
    }
  }

  # ── Summaries ──
  grid$predicted_response      <- predicted_response
  grid$ci_lower                <- apply(y_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  grid$ci_upper                <- apply(y_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  grid$predicted_concentration <- apply(x_mat, 1, stats::median, na.rm = TRUE)
  grid$se_concentration        <- apply(x_mat, 1, stats::sd, na.rm = TRUE)

  grid$pcov <- vapply(seq_len(n_grid), function(i) {
    se_i <- grid$se_concentration[i]
    if (!is.finite(se_i)) return(cv_x_max)
    raw_cv <- if (is_log_x) se_i * log(10) * 100
    else if (abs(grid$predicted_concentration[i]) > 1e-10)
      (se_i / abs(grid$predicted_concentration[i])) * 100
    else Inf
    min(raw_cv, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  grid$pcov_rmse <- vapply(seq_len(n_grid), function(i) {
    x_draws <- x_mat[i, ]
    x_true  <- grid$x_fit[i]
    ok <- is.finite(x_draws)
    if (sum(ok) < 2) return(cv_x_max)
    rmse <- sqrt(mean((x_draws[ok] - x_true)^2))
    raw_rrmse <- if (is_log_x) rmse * log(10) * 100
    else if (abs(x_true) > 1e-10) (rmse / abs(x_true)) * 100
    else Inf
    min(raw_rrmse, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  grid$pcov_pass  <- !is.na(grid$pcov) & grid$pcov < pcov_threshold
  grid$noise_mode <- if (!include_measurement_error) "curve_only"
                     else if (use_hetero) "measurement_heteroscedastic"
                     else "measurement_homoscedastic"

  grid <- curveRcore::enrich_grid_with_d2y(grid, is_log_response = is_log_response)
  grid
}


#' Back-Calculate Sample Concentrations from Posterior Draws
#'
#' For each test sample, evaluates the inverse model at every posterior draw to
#' produce a posterior distribution of predicted concentration. The variance
#' definition matches [predict_grid_bayes()] via `include_measurement_error`, so
#' the sample pcov is directly comparable to the precision profile.
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
#' @param include_measurement_error Logical. If TRUE (default) inject the SAME
#'   observation noise the grid uses before back-calculating, so samples lie on
#'   the measurement/CDAN profile. If FALSE, invert the observed response with no
#'   added noise (curve/parameter precision only). MUST match the value passed to
#'   [predict_grid_bayes()] for the same fit.
#'
#' @return Data frame with original sample columns plus prediction columns.
#'
#' @export
predict_samples_bayes <- function(samples, bayes_fit, curve_idx = 1L,
                                  response_variable,
                                  is_log_response = TRUE,
                                  n_draws = NULL, cv_x_max = 150,
                                  pcov_threshold = 20,
                                  is_log_x = TRUE,
                                  include_measurement_error = TRUE) {

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

  # Noise draws — pulled the same way the grid does, so the injected noise is
  # identical in definition. Only used when include_measurement_error = TRUE.
  sigma_obs_draws <- as.numeric(draws[["sigma_obs"]])
  nu_draws        <- as.numeric(draws[["nu"]])
  has_hetero      <- all(c("log_sigma0", "log_sigma_slope") %in% names(draws))
  log_sigma0_draws      <- if (has_hetero) as.numeric(draws[["log_sigma0"]])      else NULL
  log_sigma_slope_draws <- if (has_hetero) as.numeric(draws[["log_sigma_slope"]]) else NULL
  use_hetero <- isTRUE(bayes_fit$stan_data$use_heteroscedastic_noise == 1L)

  S <- length(a_draws)
  if (!is.null(n_draws) && n_draws < S) {
    idx <- sample.int(S, n_draws)
    a_draws <- a_draws[idx]; b_draws <- b_draws[idx]
    c_draws <- c_draws[idx]; d_draws <- d_draws[idx]
    if (!is.null(g_draws)) g_draws <- g_draws[idx]
    sigma_obs_draws <- sigma_obs_draws[idx]
    nu_draws        <- nu_draws[idx]
    if (!is.null(log_sigma0_draws)) {
      log_sigma0_draws      <- log_sigma0_draws[idx]
      log_sigma_slope_draws <- log_sigma_slope_draws[idx]
    }
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

  # Posterior matrix: rows = samples, cols = draws.
  # measurement mode: perturb the observed response with the SAME noise model the
  #   grid uses, so the sample pcov includes the measurement-error term and the
  #   points land on the profile. curve-only mode: invert the fixed observed
  #   response (parameter uncertainty only).
  x_mat <- matrix(NA_real_, nrow = n, ncol = S)
  for (s in seq_len(S)) {
    for (i in seq_len(n)) {
      y_use <- if (isTRUE(include_measurement_error)) {
        sigma_s <- .obs_noise_sigma(fit_response[i], s, use_hetero,
                                    sigma_obs_draws, log_sigma0_draws, log_sigma_slope_draws)
        fit_response[i] + sigma_s * stats::rt(1, df = nu_draws[s])
      } else {
        fit_response[i]
      }
      x_mat[i, s] <- suppressWarnings(tryCatch(inv(y_use, s), error = function(e) NA_real_))
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

  samples$pcov <- vapply(seq_len(n), function(i) {
    if (!is.finite(x_se[i])) return(cv_x_max)
    raw_cv <- if (is_log_x) x_se[i] * log(10) * 100
    else if (abs(x_median[i]) > 1e-10) (x_se[i] / abs(x_median[i])) * 100
    else Inf
    min(raw_cv, cv_x_max, na.rm = TRUE)
  }, numeric(1))

  samples$pcov_rmse <- samples$pcov
  samples$pcov_pass <- !is.na(samples$pcov) & samples$pcov < pcov_threshold
  samples$noise_mode <- if (!include_measurement_error) "curve_only"
                        else if (use_hetero) "measurement_heteroscedastic"
                        else "measurement_homoscedastic"
  samples
}
