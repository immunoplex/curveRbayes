# =============================================================================
# stan_data.R — Build the Stan data list from preprocessed stacked data
# =============================================================================


#' Build Stan Data List for a Model Family
#'
#' Converts a preprocessed stacked data frame into the named list format
#' expected by the Stan model's `data {}` block.
#'
#' @param standards Data frame. Preprocessed stacked standards with
#'   `curve_id`, `concentration`, and response column.
#' @param response_variable Character. Response column name.
#' @param priors Named list from [curveRbayes::compute_dynamic_priors()].
#' @param model_family Character. One of the curveRcore model names.
#' @param curve_id_map Named integer vector mapping curve_id values to
#'   1-based Stan indices.
#' @param blanks Data frame or NULL. Blank well data. Must contain `curve_id`
#'   and the response column (on the fitting scale). Passed directly to Stan
#'   to anchor the lower asymptote.
#' @param use_heteroscedastic_noise Logical. If TRUE, passes
#'   `use_heteroscedastic_noise = 1L` to Stan, enabling the power-of-mean
#'   noise model `sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))`.
#'   If FALSE (default), the homoscedastic `sigma_obs` path is used.
#' @param grainsize Integer. reduce_sum grain size. Default 1.
#'   NOTE: grainsize is declared in the Stan data block but reduce_sum is
#'   not yet active in the model bodies; this field is reserved for a future
#'   threading implementation and is currently ignored by Stan.
#'
#' @return Named list suitable for `cmdstanr::sample(data = ...)`.
#' @export
build_stan_data <- function(standards,
                            response_variable,
                            priors,
                            model_family = "logistic4",
                            curve_id_map,
                            blanks = NULL,
                            use_heteroscedastic_noise = FALSE,
                            grainsize = 1L) {

  n_curves <- length(curve_id_map)

  all_x     <- standards$concentration
  all_y     <- standards[[response_variable]]
  all_idx   <- as.integer(curve_id_map[as.character(standards$curve_id)])

  N_obs <- length(all_y)

  # ── Blank data ──
  if (!is.null(blanks) && nrow(blanks) > 0 &&
      response_variable %in% names(blanks) &&
      "curve_id" %in% names(blanks)) {
    blank_idx      <- as.integer(curve_id_map[as.character(blanks$curve_id)])
    blank_response <- blanks[[response_variable]]
    # Drop any blanks whose curve_id is not in the map (e.g., filtered plates)
    valid          <- !is.na(blank_idx)
    blank_idx      <- blank_idx[valid]
    blank_response <- blank_response[valid]
    N_blanks       <- length(blank_response)
  } else {
    N_blanks       <- 0L
    blank_idx      <- integer(0)
    blank_response <- numeric(0)
  }

  stan_data <- list(
    N_obs       = N_obs,
    N_plates    = n_curves,
    grainsize   = as.integer(grainsize),
    plate_idx   = all_idx,
    x           = all_x,
    y           = all_y,

    # Priors
    prior_a_mu        = priors$prior_a_mu,
    prior_a_sigma     = priors$prior_a_sigma,
    prior_d_mu        = priors$prior_d_mu,
    prior_d_sigma     = priors$prior_d_sigma,
    prior_log_b_mu    = priors$prior_log_b_mu,
    prior_log_b_sigma = priors$prior_log_b_sigma,

    # Blank data (may be empty)
    N_blanks        = N_blanks,
    blank_plate_idx = blank_idx,
    blank_response  = blank_response,

    # Noise model switch
    use_heteroscedastic_noise = as.integer(use_heteroscedastic_noise),

    # Heteroscedastic noise priors (always passed; Stan uses them only
    # when use_heteroscedastic_noise = 1)
    prior_log_sigma0_mu         = priors$prior_log_sigma0_mu,
    prior_log_sigma0_sigma      = priors$prior_log_sigma0_sigma,
    prior_log_sigma_slope_mu    = priors$prior_log_sigma_slope_mu,
    prior_log_sigma_slope_sigma = priors$prior_log_sigma_slope_sigma
  )

  # c-parameter priors differ by family
  if (model_family == "loglogistic4") {
    stan_data$prior_log_c_mu    <- priors$prior_log_c_mu
    stan_data$prior_log_c_sigma <- priors$prior_log_c_sigma
  } else {
    stan_data$prior_c_mu    <- priors$prior_c_mu
    stan_data$prior_c_sigma <- priors$prior_c_sigma
  }

  # g priors for 5-parameter families
  if (model_family %in% c("logistic5", "loglogistic5")) {
    stan_data$prior_log_g_sd       <- priors$prior_log_g_sd
    stan_data$prior_log_g_plate_sd <- priors$prior_log_g_plate_sd
  }

  stan_data
}
