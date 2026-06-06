# =============================================================================
# test-elisa-bayes.R — Bayesian tests using the ELISA example dataset
#
# ELISA data: OD response (0–4 range), 6 curve_ids, single analyte.
# Plates 1-3 generated from 5PL, plates 4-6 from Gompertz.
# Tests the hierarchical model on a different response scale than the
# bead assay, exercising medium/low scale class priors.
#
# Stan tests skip gracefully when cmdstanr is unavailable.
# =============================================================================

library(curveRcore)


# =============================================================================
# Skip helper
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
# Setup: preprocess ELISA data
# =============================================================================

data("elisa_assay_example", package = "curveRcore")

elisa_std_raw <- elisa_assay_example$standards
elisa_samples <- elisa_assay_example$samples

elisa_prepped <- curveRcore::preprocess_standards(
  data                 = elisa_std_raw,
  antigen_settings     = list(standard_curve_concentration = 10000),
  response_variable    = "od",
  independent_variable = "concentration",
  is_log_response      = TRUE,
  is_log_independent   = TRUE,
  apply_prozone        = TRUE
)

elisa_standards <- elisa_prepped$data


# =============================================================================
# Priors: ELISA scale
# =============================================================================

test_that("ELISA priors have reasonable scale for OD data", {
  priors <- compute_dynamic_priors(elisa_standards, "od")

  # OD data on log10 scale: y_min ~ log10(0.01) = -2, y_max ~ log10(3) = 0.48
  expect_true(priors$prior_a_mu < 1,
              info = paste("prior_a_mu =", priors$prior_a_mu))
  expect_true(priors$prior_d_mu > -1,
              info = paste("prior_d_mu =", priors$prior_d_mu))
  expect_true(priors$prior_a_sigma > 0)
  expect_true(priors$prior_c_sigma > 0)
})


# =============================================================================
# Stan data: ELISA
# =============================================================================

test_that("ELISA Stan data builds correctly for subset of plates", {
  # Use plates 1-3 only
  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(1, 2, 3), ]
  priors <- compute_dynamic_priors(std_sub, "od")
  curve_ids <- sort(unique(std_sub$curve_id))
  curve_id_map <- stats::setNames(seq_along(curve_ids), as.character(curve_ids))

  sdata <- build_stan_data(
    standards = std_sub,
    response_variable = "od",
    priors = priors,
    model_family = "logistic4",
    curve_id_map = curve_id_map
  )

  expect_equal(sdata$N_plates, 3)
  expect_equal(sdata$N_obs, nrow(std_sub))
  expect_true(all(sdata$plate_idx %in% 1:3))
  expect_equal(sdata$N_blanks, 0L)
})


# =============================================================================
# Bayesian fit: ELISA plates 1-3 (5PL-generated)
# =============================================================================

test_that("ELISA Bayesian fit works on plates 1-3 with logistic4", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(1, 2, 3), ]
  samp_sub <- elisa_samples[elisa_samples$curve_id %in% c(1, 2, 3), ]

  result <- fit_calibration_bayes(
    standards      = std_sub,
    samples        = samp_sub,
    response_var   = "od",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 123,
    verbose        = FALSE
  )

  # ELISA Bayesian fit works on plates 1-3 with logistic4
  # CHANGE: expect multiplate
  expect_s3_class(result, "calibration_result_multiplate")
  expect_equal(length(result$plates), 3)

  # Per-plate checks
  cr1 <- result$plates[["1"]]
  expect_equal(cr1$meta$method, "bayesian")

  tbl <- summary_table_bayes(result)  # unchanged — handles multiplate
  expect_equal(nrow(tbl), 3)

  # a checks use tbl (already per-curve), so these are fine as-is
  expect_true(all(tbl$a_mean > -3))
  expect_true(all(tbl$a_mean < 1))
  expect_true(all(tbl$d_mean > -1))
})


# =============================================================================
# Bayesian fit: ELISA plates 4-6 (Gompertz-generated)
# =============================================================================

test_that("ELISA Bayesian fit works on Gompertz-generated plates", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(4, 5, 6), ]

  result <- fit_calibration_bayes(
    standards      = std_sub,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 456,
    run_loo        = TRUE,
    verbose        = FALSE
  )

  # ELISA Bayesian fit works on Gompertz-generated plates
  expect_s3_class(result, "calibration_result_multiplate")
  cr4 <- result$plates[["4"]]
  expect_equal(length(cr4$ensemble), 2)
  expect_true(cr4$selection$best_model_name %in% c("logistic4", "gompertz4"))
})


# =============================================================================
# Bayesian fit: all 6 ELISA plates hierarchically
# =============================================================================

test_that("ELISA Bayesian fit works on all 6 plates", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  result <- fit_calibration_bayes(
    standards      = elisa_standards,
    samples        = elisa_samples,
    response_var   = "od",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 789,
    verbose        = FALSE
  )

  # ELISA Bayesian fit works on all 6 plates
  expect_s3_class(result, "calibration_result_multiplate")
  expect_equal(length(result$plates), 6)

  tbl <- summary_table_bayes(result)
  expect_equal(nrow(tbl), 6)
  expect_true(all(!is.na(tbl$a_mean)))
  expect_true(all(!is.na(tbl$d_mean)))
})


# =============================================================================
# Sample predictions: ELISA Bayesian
# =============================================================================

test_that("ELISA Bayesian sample predictions are produced", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(1, 2), ]
  samp_sub <- elisa_samples[elisa_samples$curve_id %in% c(1, 2), ]

  result <- fit_calibration_bayes(
    standards      = std_sub,
    samples        = samp_sub,
    response_var   = "od",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 101,
    verbose        = FALSE
  )

  # ELISA Bayesian sample predictions
  expect_s3_class(result, "calibration_result_multiplate")
  cr1 <- result$plates[["1"]]
  expect_true(!is.null(cr1$samples))
  expect_true(nrow(cr1$samples) > 0)
  expect_true("predicted_log10_concentration" %in% names(cr1$samples))
  expect_true("curve_id" %in% names(cr1$samples))
})


# =============================================================================
# Grid monotonicity: ELISA Bayesian
# =============================================================================

test_that("ELISA Bayesian grid predicted_response is monotonic", {
  skip_if_no_cmdstanr()
  skip_on_cran()

  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(1, 2, 3), ]

  result <- fit_calibration_bayes(
    standards      = std_sub,
    response_var   = "od",
    model_names    = "logistic4",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 200,
    sampling       = 200,
    seed           = 202,
    verbose        = FALSE
  )

  # Grid monotonicity — check per plate
  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]
    y <- cr$grid$predicted_response
    if (all(is.finite(y))) {
      expect_true(all(diff(y) >= -0.01),
                  info = paste("curve_id", cid, "not monotonic"))
    }
  }

})


# =============================================================================
# Freq vs Bayes comparison on ELISA: parameters should be close
# =============================================================================

test_that("ELISA freq and Bayes parameter estimates are comparable", {
  skip_if_not_installed("curveRfreq")
  skip_if_no_cmdstanr()
  skip_on_cran()
  library(curveRfreq)

  std_sub <- elisa_standards[elisa_standards$curve_id %in% c(1, 2, 3), ]
  samp_sub <- elisa_samples[elisa_samples$curve_id %in% c(1, 2, 3), ]

  # Frequentist
  freq_result <- curveRfreq::fit_calibration_freq_multiplate(
    standards      = std_sub,
    samples        = samp_sub,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    verbose        = FALSE
  )

  freq_tbl <- curveRfreq::summary_table(freq_result)

  # Bayesian
  bayes_result <- fit_calibration_bayes(
    standards      = std_sub,
    samples        = samp_sub,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    chains         = 2,
    warmup         = 300,
    sampling       = 300,
    seed           = 303,
    verbose        = FALSE
  )

  bayes_tbl <- summary_table_bayes(bayes_result)

  # Freq vs Bayes comparison — access Bayes per-plate
  # bayes_result is now multiplate, so access per-curve:
  # d parameter should be close between freq and Bayes
  # (d is the best-constrained parameter — both methods should agree)
  for (cid_str in as.character(1:3)) {
    # freq_d  <- freq_tbl$d[freq_tbl$curve_id == cid_str]
    # bayes_d <- bayes_tbl$d_mean[bayes_tbl$curve_id == cid_str]

    freq_d  <- freq_tbl$d[freq_tbl$curve_id == cid_str]
    bayes_d <- bayes_tbl$d_mean[bayes_tbl$curve_id == cid_str]

    if (length(freq_d) == 1 && length(bayes_d) == 1 &&
        is.finite(freq_d) && is.finite(bayes_d)) {
      expect_equal(freq_d, bayes_d, tolerance = 0.3,
                   info = paste("curve_id", cid_str,
                                "freq_d =", round(freq_d, 3),
                                "bayes_d =", round(bayes_d, 3)))
    }
  }
})
