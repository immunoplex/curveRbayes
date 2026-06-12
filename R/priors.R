# =============================================================================
# priors.R — Data-Adaptive Prior Computation
#
# Computes hyperprior parameters for the hierarchical Stan models from
# preprocessed standard curve data (already on the fitting scale).
# No blank data — priors come from the standards only.
#
# When fixed_a is supplied, the prior on a is tightened to act as a
# Bayesian soft constraint.
#
# Heteroscedastic noise priors (for use_heteroscedastic_noise = 1):
#   log_sigma0     — intercept of log(sigma) vs log(|mu|) line.
#                    Initialised from the observed residual SD.
#   log_sigma_slope — slope. Prior centred on 1 (proportional noise)
#                    with SD 0.5 to allow sub- or super-proportional noise.
# =============================================================================


#' Compute Data-Adaptive Priors for Stan Models
#'
#' Inspects preprocessed standard curve data and produces all hyperprior
#' locations and scales needed by the Stan model `data {}` block.
#'
#' @param data Data frame of preprocessed standards with `concentration`
#'   and response columns (already on the fitting scale).
#' @param response_variable Character. Response column name.
#' @param fixed_a Numeric or NULL. If non-NULL, the prior on `a` is
#'   tightened to a soft constraint centered on this value (on the
#'   fitting scale).
#' @param model_family Character. One of the curveRcore model names.
#'   Default `"logistic4"`.
#'
#' @return Named list of hyperprior values matching the Stan `data` block,
#'   including `prior_log_sigma0_mu`, `prior_log_sigma0_sigma`,
#'   `prior_log_sigma_slope_mu`, and `prior_log_sigma_slope_sigma` for
#'   the heteroscedastic noise path.
#' @export
compute_dynamic_priors <- function(data, response_variable,
                                   fixed_a = NULL,
                                   model_family = "logistic4") {

  y <- data[[response_variable]]
  y <- y[is.finite(y)]
  x <- data$concentration[is.finite(data$concentration)]

  y_min   <- min(y)
  y_max   <- max(y)
  y_range <- y_max - y_min
  x_mid   <- mean(range(x))
  x_range <- diff(range(x))

  # ── Lower asymptote (a) ──
  if (!is.null(fixed_a)) {
    # Soft constraint: tight prior centered on fixed_a
    prior_a_mu    <- fixed_a
    prior_a_sigma <- max(abs(y_range) * 0.01, 1e-4)
  } else {
    # Data-adaptive: centered on observed minimum
    prior_a_mu    <- y_min
    prior_a_sigma <- y_range * 0.3
  }

  # ── Upper asymptote (d) ──
  prior_d_mu    <- y_max + y_range * 0.1
  prior_d_sigma <- y_range * 0.3

  # ── Hill slope (b) ──
  # b > 0, estimated on log scale; median ~ 1.0
  prior_log_b_mu    <- 0.0
  prior_log_b_sigma <- 0.7

  # ── Inflection point (c) ──
  prior_c_mu    <- x_mid
  prior_c_sigma <- x_range * 0.5

  # ── Asymmetry (g) — 5-parameter models only ──
  # g = 1 is the 4PL; regularise toward 4PL
  prior_log_g_sd       <- 0.5
  prior_log_g_plate_sd <- 0.3

  # ── Heteroscedastic noise priors ──
  # log_sigma0: intercept of log(sigma) ~ log(|mu|) line.
  # Anchored to log of a robust estimate of the residual SD.
  # We use the IQR / 1.35 of y as a scale-invariant noise estimate.
  y_noise_est <- max(stats::IQR(y) / 1.35, 1e-6)
  prior_log_sigma0_mu    <- log(y_noise_est * 0.3)  # expect noise < IQR
  prior_log_sigma0_sigma <- 1.5                       # wide — data will dominate

  # log_sigma_slope: 1 = proportional noise (CV constant), 0 = additive,
  # 2 = strongly heteroscedastic. Prior centred on 1, SD = 0.5.
  prior_log_sigma_slope_mu    <- 1.0
  prior_log_sigma_slope_sigma <- 0.5

  priors <- list(
    prior_a_mu         = prior_a_mu,
    prior_a_sigma      = prior_a_sigma,
    prior_d_mu         = prior_d_mu,
    prior_d_sigma      = prior_d_sigma,
    prior_log_b_mu     = prior_log_b_mu,
    prior_log_b_sigma  = prior_log_b_sigma,
    # Heteroscedastic noise hypers (always included; used by Stan
    # when use_heteroscedastic_noise = 1)
    prior_log_sigma0_mu         = prior_log_sigma0_mu,
    prior_log_sigma0_sigma      = prior_log_sigma0_sigma,
    prior_log_sigma_slope_mu    = prior_log_sigma_slope_mu,
    prior_log_sigma_slope_sigma = prior_log_sigma_slope_sigma
  )

  # c-parameter priors differ by family
  if (model_family == "loglogistic4") {
    conc_vals <- x[is.finite(x) & x > 0]
    geo_mid <- exp(mean(log(conc_vals)))
    priors$prior_log_c_mu    <- log(geo_mid)
    priors$prior_log_c_sigma <- 1.5
  } else {
    priors$prior_c_mu    <- x_mid
    priors$prior_c_sigma <- x_range * 0.5
  }

  if (model_family %in% c("logistic5", "loglogistic5")) {
    priors$prior_log_g_sd       <- prior_log_g_sd
    priors$prior_log_g_plate_sd <- prior_log_g_plate_sd
  }

  priors
}
