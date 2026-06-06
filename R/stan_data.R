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
#' @param grainsize Integer. reduce_sum grain size. Default 1.
#'
#' @return Named list suitable for `cmdstanr::sample(data = ...)`.
#' @export
build_stan_data <- function(standards,
                            response_variable,
                            priors,
                            model_family = "logistic4",
                            curve_id_map,
                            grainsize = 1L) {

  n_curves <- length(curve_id_map)

  all_x     <- standards$concentration
  all_y     <- standards[[response_variable]]
  all_idx   <- as.integer(curve_id_map[as.character(standards$curve_id)])

  N_obs <- length(all_y)

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

    # No blanks
    N_blanks        = 0L,
    blank_plate_idx = integer(0),
    blank_response  = numeric(0)
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
