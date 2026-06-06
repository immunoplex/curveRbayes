# =============================================================================
# loo_selection.R — PSIS-LOO cross-validation and Bayesian stacking weights
# =============================================================================

#' Compute LOO-CV for a Fitted Bayesian Model
#'
#' Extracts the log_lik generated quantity and computes PSIS-LOO.
#'
#' @param bayes_fit Output of [curveRbayes::fit_bayes_single()].
#' @return A `loo` object (from the loo package).
#' @export
compute_loo <- function(bayes_fit) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("loo package is required for model comparison")

  draws <- bayes_fit$draws
  ll_cols <- grep("^log_lik\\[", names(draws), value = TRUE)
  if (length(ll_cols) == 0)
    stop("No log_lik columns found in posterior draws")

  # Convert to plain matrix to avoid posterior::draws_df subsetting warnings
  log_lik_matrix <- as.matrix(draws)[, ll_cols, drop = FALSE]
  suppressWarnings(loo::loo(log_lik_matrix))
}


#' Compare Models via LOO-CV and Stacking Weights
#'
#' Given a list of fitted Bayesian models (one per model family),
#' computes LOO for each and Bayesian stacking weights.
#'
#' @param fits Named list of [curveRbayes::fit_bayes_single()] outputs.
#' @return Named list with `best_model_name`, `criterion`, `loo_results`,
#'   `comparison`, `weights`.
#' @export
compare_models_loo <- function(fits) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("loo package is required for model comparison")

  loo_list <- lapply(fits, compute_loo)
  names(loo_list) <- names(fits)

  comp <- loo::loo_compare(loo_list)

  # Stacking weights
  stacking_wts <- tryCatch(
    suppressWarnings(loo::stacking_weights(loo_list)),
    error = function(e) {
      w <- rep(0, length(fits))
      names(w) <- names(fits)
      w[rownames(comp)[1]] <- 1
      w
    }
  )

  best_name <- rownames(comp)[1]

  list(
    best_model_name = best_name,
    criterion       = "LOO",
    loo_results     = loo_list,
    comparison      = comp,
    weights         = as.numeric(stacking_wts)
  )
}
