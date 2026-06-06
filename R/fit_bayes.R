# =============================================================================
# fit_bayes.R — Core Bayesian fitting wrapper
#
# Compiles the Stan model (once, cached), runs MCMC via cmdstanr,
# and packages the result for downstream use.
# =============================================================================


#' Get Path to a curveRbayes Stan Model File
#'
#' @param model_family Character. One of the curveRcore model names.
#' @return Absolute path to the `.stan` file.
#' @keywords internal
stan_model_path <- function(model_family) {
  fname <- switch(model_family,
                  logistic4    = "hierarchical_logistic4.stan",
                  logistic5    = "hierarchical_logistic5.stan",
                  loglogistic4 = "hierarchical_loglogistic4.stan",
                  loglogistic5 = "hierarchical_loglogistic5.stan",
                  gompertz4    = "hierarchical_gompertz4.stan",
                  stop("Unknown model_family: ", model_family)
  )
  system.file("stan", fname, package = "curveRbayes", mustWork = TRUE)
}


#' Compile a curveRbayes Stan Model (Cached)
#'
#' Compiles the Stan model via cmdstanr. Compilation is cached by
#' cmdstanr so subsequent calls are instant.
#'
#' @param model_family Character. One of the curveRcore model names.
#' @return A `CmdStanModel` object.
#' @export
compile_stan_model <- function(model_family = "logistic4") {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("cmdstanr is required for curveRbayes. Install via:\n",
         "  install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))")

  path <- stan_model_path(model_family)
  cmdstanr::cmdstan_model(path)
}


#' Fit a Single Model Family via MCMC
#'
#' Runs HMC/NUTS sampling for one model family across all curve_ids
#' in the data. Returns a structured list containing the CmdStanMCMC
#' fit, posterior draws, and metadata.
#'
#' @param stan_data Named list from [curveRbayes::build_stan_data()].
#' @param model_family Character. Model name.
#' @param chains Integer. Default 4.
#' @param warmup Integer. Default 1000.
#' @param sampling Integer. Default 1000.
#' @param adapt_delta Numeric. Default 0.9.
#' @param max_treedepth Integer. Default 12.
#' @param seed Integer or NULL.
#' @param compiled_model Optional pre-compiled CmdStanModel.
#' @param verbose Logical.
#'
#' @return A named list with `model_family`, `fit`, `draws`,
#'   `n_curves`, `stan_data`, `diagnostics`.
#'
#' @export
fit_bayes_single <- function(stan_data,
                             model_family = "logistic4",
                             chains = 4L,
                             warmup = 1000L,
                             sampling = 1000L,
                             adapt_delta = 0.9,
                             max_treedepth = 12L,
                             seed = NULL,
                             compiled_model = NULL,
                             verbose = FALSE) {

  mod <- compiled_model %||% compile_stan_model(model_family)

  if (verbose) message("[fit_bayes] Sampling ", model_family,
                       " (", chains, " chains \u00d7 ", sampling, " draws) ...")

  mcmc_fit <- mod$sample(
    data            = stan_data,
    chains          = chains,
    iter_warmup     = warmup,
    iter_sampling   = sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    seed            = seed,
    refresh         = if (verbose) 200 else 0,
    show_messages   = verbose,
    show_exceptions = verbose
  )

  draws <- posterior::as_draws_df(mcmc_fit$draws())

  diag <- list(
    num_divergent     = sum(mcmc_fit$diagnostic_summary()$num_divergent),
    num_max_treedepth = sum(mcmc_fit$diagnostic_summary()$num_max_treedepth),
    ebfmi             = mcmc_fit$diagnostic_summary()$ebfmi
  )

  if (verbose) {
    message("[fit_bayes] Done. Divergences: ", diag$num_divergent,
            "  Max treedepth: ", diag$num_max_treedepth)
  }

  list(
    model_family = model_family,
    fit          = mcmc_fit,
    draws        = draws,
    n_curves     = stan_data$N_plates,
    stan_data    = stan_data,
    diagnostics  = diag
  )
}


#' Extract Curve-Level Posterior Summaries
#'
#' Computes posterior mean, SD, and quantiles for one curve's parameters.
#'
#' @param bayes_fit Output of [curveRbayes::fit_bayes_single()].
#' @param curve_idx Integer. Which curve (1-based Stan index).
#' @param probs Numeric vector of quantiles. Default c(0.025, 0.5, 0.975).
#'
#' @return Data frame with columns: term, mean, sd, q2.5, q50, q97.5.
#' @export
extract_curve_params <- function(bayes_fit, curve_idx = 1L,
                                 probs = c(0.025, 0.5, 0.975)) {

  draws  <- bayes_fit$draws
  family <- bayes_fit$model_family
  p      <- curve_idx

  param_names <- if (family %in% c("logistic5", "loglogistic5")) {
    c("a", "b", "c_par", "d", "g")
  } else {
    c("a", "b", "c_par", "d")
  }

  stan_names <- paste0(param_names, "[", p, "]")
  present <- stan_names %in% names(draws)

  if (!any(present)) {
    warning("No parameters found for curve index ", p)
    return(data.frame(term = character(), mean = numeric(),
                      sd = numeric(), stringsAsFactors = FALSE))
  }

  rows <- lapply(seq_along(param_names), function(j) {
    sn <- stan_names[j]
    if (!(sn %in% names(draws))) return(NULL)
    vals <- as.numeric(draws[[sn]])
    qs <- stats::quantile(vals, probs = probs)
    data.frame(
      term  = sub("c_par", "c", param_names[j]),
      mean  = mean(vals),
      sd    = stats::sd(vals),
      q2.5  = qs[1],
      q50   = qs[2],
      q97.5 = qs[3],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}
