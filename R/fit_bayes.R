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


# =============================================================================
# Schema-expansion extractors (calib_hyperparam / calib_draws / calib_fit_diag)
#
# These read the `bf` object returned by fit_bayes_single() — specifically
# bf$draws (posterior::as_draws_df, all params incl. transformed b/g and the
# non-indexed population/noise scalars) and bf$fit (the CmdStanMCMC object).
# They fill the curveRcore `population` slot contract:
#   params   : data.frame(term, estimate, std_error, q_lo, q_med, q_hi)
#   draws    : named list term -> numeric() (iteration-ordered; shared order)
#   fit_diag : named list of per-fit sampler diagnostics
# =============================================================================

# Population/noise scalar term names in the draws = columns that are NOT
# plate-indexed ("[..]"), NOT non-centered auxiliaries (raw_*), and NOT
# posterior bookkeeping (.chain/.iteration/.draw, lp__). Auto-discovers the
# exact .stan parameters{} scalars (mu_*, sigma_*, sigma_obs, nu, sigma_blank,
# log_sigma0, log_sigma_slope) for whichever model was fitted.
.population_terms <- function(draws) {
  nm <- names(draws)
  meta    <- nm %in% c(".chain", ".iteration", ".draw", ".draw_id") | nm == "lp__"
  indexed <- grepl("\\[", nm)
  raw     <- grepl("^raw_", nm)
  nm[!meta & !indexed & !raw]
}

#' Summarise population / noise parameters from a Bayesian fit
#'
#' @param bayes_fit Output of [fit_bayes_single()].
#' @param probs Quantiles. Default c(0.025, 0.5, 0.975).
#' @return data.frame(term, estimate, std_error, q_lo, q_med, q_hi) — the
#'   curveRcore `population$params` shape. One row per group-level scalar.
#' @export
extract_population_summary <- function(bayes_fit, probs = c(0.025, 0.5, 0.975)) {
  draws <- bayes_fit$draws
  terms <- .population_terms(draws)
  if (length(terms) == 0L) return(data.frame())
  rows <- lapply(terms, function(t) {
    v  <- as.numeric(draws[[t]])
    qs <- stats::quantile(v, probs = probs, names = FALSE, na.rm = TRUE)
    data.frame(term = t, estimate = mean(v), std_error = stats::sd(v),
               q_lo = qs[1], q_med = qs[2], q_hi = qs[3],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Posterior draw vectors for the population / noise scalars
#'
#' @param bayes_fit Output of [fit_bayes_single()].
#' @return Named list `term -> numeric()`, iteration-ordered (all terms share
#'   one order, so they can be column-bound into a joint posterior). For
#'   calib_draws (param_scope = "population"); emit only when persist_draws.
#' @export
extract_population_draws <- function(bayes_fit) {
  draws <- bayes_fit$draws
  terms <- .population_terms(draws)
  stats::setNames(lapply(terms, function(t) as.numeric(draws[[t]])), terms)
}

#' Posterior draw vectors for one plate's curve parameters (a,b,c,d,g)
#'
#' @param bayes_fit Output of [fit_bayes_single()].
#' @param curve_idx Integer plate index (1-based Stan index).
#' @return Named list `term -> numeric()` (c_par emitted as c). For calib_draws
#'   (param_scope = "curve"); emit only when persist_draws.
#' @export
extract_curve_draws <- function(bayes_fit, curve_idx) {
  draws  <- bayes_fit$draws
  family <- bayes_fit$model_family
  pnames <- if (family %in% c("logistic5", "loglogistic5"))
    c("a", "b", "c_par", "d", "g") else c("a", "b", "c_par", "d")
  out <- list()
  for (pn in pnames) {
    sn <- paste0(pn, "[", curve_idx, "]")
    if (sn %in% names(draws)) out[[sub("c_par", "c", pn)]] <- as.numeric(draws[[sn]])
  }
  out
}

#' Per-fit sampler diagnostics from a Bayesian fit
#'
#' @param bayes_fit Output of [fit_bayes_single()].
#' @return Named list matching the calib_fit_diag Bayesian columns. `converged`
#'   is left NULL here (the caller supplies it from the ensemble's converged
#'   flag); rhat/ess come from the CmdStanMCMC $summary(), divergences/treedepth/
#'   ebfmi from bf$diagnostics, timing from $time(), seed from $metadata().
#' @export
extract_fit_diag <- function(bayes_fit) {
  fit  <- bayes_fit$fit
  diag <- bayes_fit$diagnostics %||% list()

  rhat_max <- ess_bulk_min <- ess_tail_min <- NA_real_
  smry <- tryCatch(fit$summary(), error = function(e) NULL)
  if (!is.null(smry)) {
    if ("rhat"     %in% names(smry)) rhat_max     <- suppressWarnings(max(smry$rhat,     na.rm = TRUE))
    if ("ess_bulk" %in% names(smry)) ess_bulk_min <- suppressWarnings(min(smry$ess_bulk, na.rm = TRUE))
    if ("ess_tail" %in% names(smry)) ess_tail_min <- suppressWarnings(min(smry$ess_tail, na.rm = TRUE))
  }
  fin <- function(x) if (length(x) && is.finite(x)) x else NA_real_

  n_draws_total <- tryCatch(nrow(bayes_fit$draws), error = function(e) NA_integer_)
  n_div <- diag$num_divergent %||% NA_integer_
  fit_seconds <- tryCatch({ t <- fit$time(); t$total %||% NA_real_ }, error = function(e) NA_real_)
  fit_seed <- tryCatch({ s <- fit$metadata()$seed; if (length(s)) as.numeric(s)[1] else NA_real_ },
                       error = function(e) NA_real_)
  ebfmi_min <- tryCatch(suppressWarnings(min(diag$ebfmi, na.rm = TRUE)), error = function(e) NA_real_)

  list(
    fit_seconds       = fin(fit_seconds),
    n_iterations      = NA_integer_,               # MCMC: warmup+sampling live in meta
    converged         = NULL,                       # set by caller (ensemble flag)
    fit_seed          = if (is.finite(fit_seed)) fit_seed else NA_real_,
    rhat_max          = fin(rhat_max),
    ess_bulk_min      = fin(ess_bulk_min),
    ess_tail_min      = fin(ess_tail_min),
    n_divergent       = n_div,
    pct_divergent     = if (!is.na(n_div) && !is.na(n_draws_total) && n_draws_total > 0)
                          100 * n_div / n_draws_total else NA_real_,
    max_treedepth_hit = diag$num_max_treedepth %||% NA_integer_,
    ebfmi_min         = fin(ebfmi_min)
  )
}
