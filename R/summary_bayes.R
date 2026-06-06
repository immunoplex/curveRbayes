# =============================================================================
# summary_bayes.R — Summary extraction for Bayesian results
# =============================================================================


#' Extract a Per-Curve Summary Table from Bayesian Results
#'
#' Works with both `calibration_result_multiplate` (new format) and
#' single `calibration_result` (legacy format).
#'
#' @param cr A `calibration_result_multiplate` or `calibration_result`
#'   with `method = "bayesian"`.
#'
#' @return Data frame with one row per curve_id.
#'
#' @export
summary_table_bayes <- function(cr) {

  # Handle multiplate format (new)
  if (inherits(cr, "calibration_result_multiplate")) {
    rows <- lapply(seq_along(cr$plates), function(i) {
      plate <- cr$plates[[i]]
      cid   <- names(cr$plates)[i]
      if (is.null(plate)) {
        return(data.frame(curve_id = cid, best_model = NA_character_,
                          stringsAsFactors = FALSE))
      }

      best <- plate$selection$best_model_name
      ens  <- if (!is.na(best)) plate$ensemble[[best]] else NULL
      pf   <- if (!is.null(ens)) ens$parameters else NULL

      row <- data.frame(curve_id = cid, best_model = best,
                        stringsAsFactors = FALSE)

      if (!is.null(pf) && is.data.frame(pf) && nrow(pf) > 0) {
        for (j in seq_len(nrow(pf))) {
          term <- pf$term[j]
          row[[paste0(term, "_mean")]] <- pf$mean[j]
          row[[paste0(term, "_sd")]]   <- pf$sd[j]
        }
      }

      if (!is.null(ens$fit_stats)) {
        row$n_divergent     <- ens$fit_stats$n_divergent
        row$n_max_treedepth <- ens$fit_stats$n_max_treedepth
      }

      row
    })
    return(do.call(rbind, rows))
  }

  # Handle single calibration_result (legacy format)
  stopifnot(inherits(cr, "calibration_result"))

  best_name <- cr$selection$best_model_name
  if (is.na(best_name) || is.null(cr$ensemble[[best_name]])) {
    return(data.frame(curve_id = character(), best_model = character(),
                      stringsAsFactors = FALSE))
  }

  ens         <- cr$ensemble[[best_name]]
  params_list <- ens$parameters

  # If parameters is a data frame (single curve), wrap in list
  if (is.data.frame(params_list)) {
    cid <- cr$meta$curve_id
    params_list <- stats::setNames(list(params_list), cid)
  }

  n_curves  <- length(params_list)
  curve_ids <- names(params_list)

  rows <- lapply(seq_len(n_curves), function(idx) {
    cid <- if (!is.null(curve_ids)) curve_ids[idx] else as.character(idx)
    pf  <- params_list[[idx]]

    row <- data.frame(curve_id = cid, best_model = best_name,
                      stringsAsFactors = FALSE)

    if (!is.null(pf) && is.data.frame(pf) && nrow(pf) > 0) {
      for (j in seq_len(nrow(pf))) {
        term <- pf$term[j]
        row[[paste0(term, "_mean")]] <- pf$mean[j]
        row[[paste0(term, "_sd")]]   <- pf$sd[j]
      }
    }

    if (!is.null(ens$fit_stats)) {
      row$n_divergent     <- ens$fit_stats$n_divergent
      row$n_max_treedepth <- ens$fit_stats$n_max_treedepth
    }

    row
  })

  do.call(rbind, rows)
}


#' Collect All Sample Predictions from Bayesian Results
#'
#' @param cr A `calibration_result_multiplate` or `calibration_result`.
#'
#' @return Data frame, or NULL if no samples.
#' @export
collect_samples_bayes <- function(cr) {

  # Handle multiplate format
  if (inherits(cr, "calibration_result_multiplate")) {
    dfs <- lapply(seq_along(cr$plates), function(i) {
      plate <- cr$plates[[i]]
      if (is.null(plate) || is.null(plate$samples) || nrow(plate$samples) == 0)
        return(NULL)
      s <- plate$samples
      s$curve_id <- names(cr$plates)[i]
      s
    })
    dfs <- Filter(Negate(is.null), dfs)
    if (length(dfs) == 0) return(NULL)
    return(do.call(rbind, dfs))
  }

  # Handle single calibration_result
  stopifnot(inherits(cr, "calibration_result"))
  if (is.null(cr$samples) || nrow(cr$samples) == 0) return(NULL)
  s <- cr$samples
  s$best_model <- cr$selection$best_model_name
  s
}
