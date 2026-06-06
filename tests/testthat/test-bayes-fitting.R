# test-bayes-fitting.R — Tests for curveRbayes fitting pipeline
#
# Uses bead_assay_example data from curveRcore, preprocessed upstream.
# Stan compilation is slow, so we use minimal chains/iterations.
# Tests marked with skip_on_cran() for CI environments without cmdstanr.

library(curveRcore)

# =============================================================================
# Skip if cmdstanr not available
# =============================================================================
skip_if_no_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    skip("cmdstanr not available")
  }
  tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) skip("CmdStan not installed")
  )
}

# =============================================================================
# Test data setup
# =============================================================================
data("bead_assay_example", package = "curveRcore")

std_raw <- bead_assay_example$standards[
  bead_assay_example$standards$curve_id %in% c(1, 2, 3), ]

prepped <- curveRcore::preprocess_standards(
  data                 = std_raw,
  antigen_settings     = list(standard_curve_concentration = 10000),
  response_variable    = "mfi",
  independent_variable = "concentration",
  is_log_response      = TRUE,
  is_log_independent   = TRUE,
  apply_prozone        = TRUE
)

standards <- prepped$data

samples <- bead_assay_example$samples[
  bead_assay_example$samples$curve_id %in% c(1, 2, 3), ]


# =============================================================================
# Priors (unchanged — these don't depend on return format)
# =============================================================================
test_that("compute_dynamic_priors returns expected fields (a free)", {
  priors <- compute_dynamic_priors(standards, "mfi")
  expect_true(is.list(priors))
  expect_true("prior_a_mu" %in% names(priors))
  expect_true("prior_d_mu" %in% names(priors))
  expect_true("prior_log_b_mu" %in% names(priors))
  expect_true("prior_c_mu" %in% names(priors))
  expect_true(priors$prior_a_sigma > 0)
})

test_that("compute_dynamic_priors with fixed_a tightens prior", {
  priors_free  <- compute_dynamic_priors(standards, "mfi", fixed_a = NULL)
  priors_fixed <- compute_dynamic_priors(standards, "mfi", fixed_a = 1.5)
  expect_equal(priors_fixed$prior_a_mu, 1.5)
  expect_true(priors_fixed$prior_a_sigma < priors_free$prior_a_sigma)
})

test_that("compute_dynamic_priors adds g priors for 5-param models", {
  priors_4 <- compute_dynamic_priors(standards, "mfi", model_family = "logistic4")
  priors_5 <- compute_dynamic_priors(standards, "mfi", model_family = "logistic5")
  expect_false("prior_log_g_sd" %in% names(priors_4))
  expect_true("prior_log_g_sd" %in% names(priors_5))
})


# =============================================================================
# Stan data construction (unchanged — operates before fitting)
# =============================================================================
test_that("build_stan_data creates valid Stan data list", {
  priors <- compute_dynamic_priors(standards, "mfi")
  curve_ids <- sort(unique(standards$curve_id))
  curve_id_map <- stats::setNames(seq_along(curve_ids), as.character(curve_ids))

  sdata <- build_stan_data(
    standards = standards,
    response_variable = "mfi",
    priors = priors,
    model_family = "logistic4",
    curve_id_map = curve_id_map
  )

  expect_equal(sdata$N_obs, nrow(standards))
  expect_equal(sdata$N_plates, 3)
  expect_equal(length(sdata$plate_idx), nrow(standards))
  expect_true(all(sdata$plate_idx %in% 1:3))
  expect_equal(sdata$N_blanks, 0L)
})


# =============================================================================
# Stan model compilation (unchanged)
# =============================================================================
test_that("compile_stan_model works for all five families", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  for (fam in curveRcore::available_models()) {
    mod <- compile_stan_model(fam)
    expect_true(inherits(mod, "CmdStanModel"),
                info = paste("Failed to compile:", fam))
  }
})


# =============================================================================
# Full Bayesian fit: logistic4 only (fast test)
# CHANGED: expect multiplate, access per-plate results
# =============================================================================
test_that("fit_calibration_bayes works with logistic4", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    samples        = samples,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    verbose        = FALSE
  )

  # CHANGED: now returns multiplate
  expect_s3_class(result, "calibration_result_multiplate")
  expect_equal(length(result$plates), 3)
  expect_true(all(c("1", "2", "3") %in% names(result$plates)))

  # Access per-plate result
  cr1 <- result$plates[["1"]]
  expect_s3_class(cr1, "calibration_result")
  expect_equal(cr1$meta$method, "bayesian")
  expect_equal(cr1$meta$curve_id, "1")
  expect_true(nrow(cr1$grid) > 0)
  expect_true(!is.null(cr1$samples))
})


# =============================================================================
# Full Bayesian fit: two models with LOO
# CHANGED: multiplate access
# =============================================================================
test_that("fit_calibration_bayes with LOO selects best model", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    run_loo        = TRUE,
    verbose        = FALSE
  )

  cr1 <- result$plates[["1"]]

  # Criterion now reflects eligibility gating
  expect_equal(cr1$selection$criterion, "LOO+eligibility")
  expect_true(cr1$selection$best_model_name %in% c("logistic4", "gompertz4"))
  expect_equal(length(cr1$ensemble), 2)

  # Eligibility fields present
  expect_true("assessments" %in% names(cr1$selection))
  expect_true("eligible_models" %in% names(cr1$selection))
  expect_false(is.null(cr1$selection$fallback))
})


# =============================================================================
# Summary table
# CHANGED: summary_table_bayes should work on multiplate directly
# =============================================================================
test_that("summary_table_bayes returns one row per curve_id", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42
  )

  tbl <- summary_table_bayes(result)
  expect_equal(nrow(tbl), 3)
  expect_true("curve_id" %in% names(tbl))
  expect_true("a_mean" %in% names(tbl))
  expect_true("b_mean" %in% names(tbl))
  # No antigen column
  expect_false("antigen" %in% names(tbl))
})


# =============================================================================
# Parameter sanity: a should be reasonable
# CHANGED: access per-plate summary
# =============================================================================
test_that("Bayesian a estimate is in reasonable range", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 500,
    sampling       = 500,
    seed           = 42
  )

  tbl <- summary_table_bayes(result)
  # a should be between 0 and 3 on log10(MFI) scale
  expect_true(all(tbl$a_mean > 0), info = paste("a_mean:", tbl$a_mean))
  expect_true(all(tbl$a_mean < 3))
  # d should be between 3 and 5
  expect_true(all(tbl$d_mean > 3))
  expect_true(all(tbl$d_mean < 5))
})


# =============================================================================
# fixed_a tightens the a posterior
# CHANGED: multiplate return
# =============================================================================
test_that("fixed_a produces tight posterior on a", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result_free <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    fixed_a        = NULL,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42
  )

  result_fixed <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    fixed_a        = 1.5,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42
  )

  tbl_free  <- summary_table_bayes(result_free)
  tbl_fixed <- summary_table_bayes(result_fixed)

  # Fixed-a posterior should have much smaller SD
  expect_true(all(tbl_fixed$a_sd < tbl_free$a_sd))
  # Fixed-a posterior mean should be near 1.5
  expect_true(all(abs(tbl_fixed$a_mean - 1.5) < 0.5))
})


# =============================================================================
# Error handling (unchanged — errors fire before return)
# =============================================================================
test_that("missing response_var raises error", {
  expect_error(
    fit_calibration_bayes(
      standards      = standards,
      response_var   = "nonexistent",
      model_names    = "logistic4",
      std_curve_conc = 10000
    ),
    "not found"
  )
})

test_that("missing curve_id raises error", {
  bad_df <- standards
  bad_df$curve_id <- NULL
  expect_error(
    fit_calibration_bayes(
      standards      = bad_df,
      response_var   = "mfi",
      model_names    = "logistic4",
      std_curve_conc = 10000
    ),
    "curve_id"
  )
})

# =============================================================================
# Per-model CDAN grids: compute_all_grids = TRUE
# =============================================================================
test_that("compute_all_grids produces per-model grids for all ensemble models", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    n_draws_predict  = 200,
    n_draws_ensemble = 100,
    compute_all_grids = TRUE,
    verbose        = FALSE
  )

  expect_s3_class(result, "calibration_result_multiplate")

  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]

    # Both models should have grids
    for (nm in names(cr$ensemble)) {
      ens <- cr$ensemble[[nm]]
      expect_true(!is.null(ens$grid),
                  info = paste("No grid for", nm, "on curve_id", cid))
      expect_true(nrow(ens$grid) > 0,
                  info = paste("Empty grid for", nm, "on curve_id", cid))
      expect_true("pcov" %in% names(ens$grid),
                  info = paste("No pcov in grid for", nm, "on curve_id", cid))
      expect_true("pcov_rmse" %in% names(ens$grid),
                  info = paste("No pcov_rmse in grid for", nm, "on curve_id", cid))
      expect_true("predicted_response" %in% names(ens$grid),
                  info = paste("No predicted_response in grid for", nm,
                               "on curve_id", cid))
    }

    # Plate-level grid should match best-model ensemble grid
    best <- cr$selection$best_model_name
    expect_equal(cr$grid$predicted_response,
                 cr$ensemble[[best]]$grid$predicted_response,
                 info = paste("Plate grid != best ensemble grid for curve_id", cid))
  }
})


test_that("per-model grids have monotonic predicted_response", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    n_draws_predict  = 200,
    n_draws_ensemble = 100,
    compute_all_grids = TRUE,
    verbose        = FALSE
  )

  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]
    for (nm in names(cr$ensemble)) {
      g <- cr$ensemble[[nm]]$grid
      if (!is.null(g)) {
        y <- g$predicted_response
        if (all(is.finite(y))) {
          expect_true(all(diff(y) >= -0.01),
                      info = paste("Non-monotonic:", nm, "curve_id", cid))
        }
      }
    }
  }
})


test_that("per-model pcov_pass is consistent with pcov < pcov_threshold", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    n_draws_predict  = 200,
    n_draws_ensemble = 100,
    compute_all_grids = TRUE,
    verbose        = FALSE
  )

  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]
    for (nm in names(cr$ensemble)) {
      g <- cr$ensemble[[nm]]$grid
      if (!is.null(g)) {
        finite_pcov <- !is.na(g$pcov)
        if (any(finite_pcov)) {
          # pcov_pass uses pcov_threshold (precision budget), not cv_x_max (hard cap)
          expected_pass <- g$pcov < 20
          expect_equal(g$pcov_pass[finite_pcov],
                       expected_pass[finite_pcov],
                       info = paste("pcov_pass wrong:", nm, "curve_id", cid))
        }
      }
    }
  }
})


test_that("compute_all_grids is forced TRUE for multi-model fits (eligibility requires it)", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  # When multiple models are fitted, compute_all_grids is forced TRUE
  # so that the eligibility gate has pcov profiles for every model.
  # The user-supplied compute_all_grids = FALSE is overridden.
  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    compute_all_grids = FALSE,   # will be overridden internally
    verbose        = FALSE
  )

  cr1 <- result$plates[["1"]]

  # All converged models have grids (forced for eligibility)
  for (nm in names(cr1$ensemble)) {
    ens <- cr1$ensemble[[nm]]
    if (!isTRUE(ens$converged)) next
    expect_true(!is.null(ens$grid),
                info = paste("Expected grid for", nm))
    expect_true("pcov" %in% names(ens$grid),
                info = paste("Expected pcov in grid for", nm))
  }

  # Plate-level grid is always the best-model grid
  expect_true(nrow(cr1$grid) > 0)
  expect_true("pcov" %in% names(cr1$grid))
})

test_that("compute_all_grids = FALSE suppresses non-best grids for single-model fit", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  # With a single model there is nothing to force — compute_all_grids = FALSE
  # has no effect but the single model still gets its grid (it IS the best).
  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    compute_all_grids = FALSE,
    verbose        = FALSE
  )

  cr1 <- result$plates[["1"]]

  # Single model: its grid exists and is both the plate grid and ensemble grid
  expect_true(nrow(cr1$grid) > 0)
  expect_true(!is.null(cr1$ensemble$logistic4$grid))
})


test_that("meta records compute_all_grids provenance", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    n_draws_predict  = 300,
    n_draws_ensemble = 150,
    compute_all_grids = TRUE,
    verbose        = FALSE
  )

  cr1 <- result$plates[["1"]]
  expect_true(cr1$meta$compute_all_grids)
  expect_equal(cr1$meta$n_draws_predict, 300)
  expect_equal(cr1$meta$n_draws_ensemble, 150)

  # Multiplate meta also records it
  expect_true(result$meta$compute_all_grids)
})


test_that("single-model fit with compute_all_grids = TRUE still works", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = standards,
    response_var   = "mfi",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 42,
    compute_all_grids = TRUE,
    verbose        = FALSE
  )

  expect_s3_class(result, "calibration_result_multiplate")

  # With one model, the best IS the only model — grid should exist
  cr1 <- result$plates[["1"]]
  expect_true(!is.null(cr1$ensemble$logistic4$grid))
  expect_true(nrow(cr1$ensemble$logistic4$grid) > 0)
})
