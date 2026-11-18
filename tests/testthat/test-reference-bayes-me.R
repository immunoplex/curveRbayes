# =============================================================================
# test-reference-bayes-me.R
#
# Regression test for the DEFAULT back-calculation: measurement noise included
# (include_measurement_error = TRUE) and seeded (0.5.2). Same fit as
# test-reference-bayes.R (curve_ids 1-3, same seed), compared with
# fixtures/reference_bayes_curve123_me.rds.
#
# Generate the reference with:  Rscript data-raw/generate_reference_bayes_me.R
#
# Tolerances. On the machine and CmdStan version that made the snapshot, the
# same seed reproduces it exactly. Elsewhere MCMC output differs like a
# different seed, so tolerances follow the run-to-run spread measured with
# different seeds (2026-11-18): inside the dynamic range estimates agree within
# a few percent; above ULOQ, where the curve is flat, they can differ by about
# 0.1 log10. The wider guard there still catches a systematic change (the
# measurement-noise change shifted those samples by up to 0.5 log10).
# =============================================================================

library(curveRcore)

TOL_PARAM_MEAN   <- 0.05
TOL_PARAM_SD     <- 0.10
TOL_RESPONSE     <- 0.02
TOL_IN_RANGE_REL <- 0.05   # final_concentration, samples in the dynamic range
TOL_OUT_LOG10    <- 0.25   # predicted_log10_concentration, samples outside it
MIN_CLASS_AGREE  <- 0.90   # share of samples with the same gate and censor class
EXPECTED_CIDS    <- c("1", "2", "3")

cmdstanr_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    tryCatch({ cmdstanr::cmdstan_version(); TRUE }, error = function(e) FALSE)
}
ref_me_path <- function() testthat::test_path("fixtures", "reference_bayes_curve123_me.rds")
ref_june_path <- function() testthat::test_path("fixtures", "reference_bayes_curve123.rds")

# The engine does not classify precision gates; curveRcore does, after the fit
# (as the worker does). Apply the same steps to the result and the snapshot.
post <- function(mp) curveRcore::classify_pcov_gates_multiplate(
  curveRcore::compute_detection_limits_multiplate(mp))

.me_cache <- new.env()
run_me <- function() {
  if (!is.null(.me_cache$result)) return(.me_cache$result)
  data("bead_assay_example", package = "curveRcore", envir = environment())
  std_raw <- bead_assay_example$standards[bead_assay_example$standards$curve_id %in% 1:3, ]
  prepped <- curveRcore::preprocess_standards(
    data = std_raw, antigen_settings = list(standard_curve_concentration = 10000),
    response_variable = "mfi", independent_variable = "concentration",
    is_log_response = TRUE, is_log_independent = TRUE, apply_prozone = TRUE)
  samples <- bead_assay_example$samples[bead_assay_example$samples$curve_id %in% 1:3, ]
  .me_cache$result <- curveRbayes::fit_calibration_bayes(
    standards = prepped$data, samples = samples, response_var = "mfi",
    model_names = c("logistic4", "gompertz4"),
    is_log_response = TRUE, is_log_independent = TRUE, std_curve_conc = 10000,
    chains = 4, warmup = 1000, sampling = 1000, adapt_delta = 0.9, seed = 20260529,
    run_loo = TRUE, n_draws_predict = 500, include_measurement_error = TRUE,
    verbose = FALSE)
  .me_cache$result
}

skip_me <- function() {
  skip_if_not(cmdstanr_available(), "cmdstanr/CmdStan not available")
  skip_if_not(file.exists(ref_me_path()),
              "Reference not found. Run: Rscript data-raw/generate_reference_bayes_me.R")
}

test_that("measurement-noise reference: same best model and parameters", {
  skip_me()
  ref <- readRDS(ref_me_path()); res <- run_me()
  for (cid in EXPECTED_CIDS) {
    cr <- res$plates[[cid]]; rr <- ref$plates[[cid]]
    expect_equal(cr$selection$best_model_name, rr$selection$best_model_name, info = cid)
    b <- cr$selection$best_model_name
    p <- cr$ensemble[[b]]$parameters; q <- rr$ensemble[[b]]$parameters
    p <- p[order(p$term), ]; q <- q[order(q$term), ]
    expect_equal(p$term, q$term, info = cid)
    expect_equal(p$mean, q$mean, tolerance = TOL_PARAM_MEAN, info = cid)
    expect_equal(p$sd, q$sd, tolerance = TOL_PARAM_SD, info = cid)
  }
})

test_that("measurement-noise reference: grid response and noise mode", {
  skip_me()
  ref <- readRDS(ref_me_path()); res <- run_me()
  for (cid in EXPECTED_CIDS) {
    g <- res$plates[[cid]]$grid; h <- ref$plates[[cid]]$grid
    expect_equal(g$predicted_response, h$predicted_response, tolerance = TOL_RESPONSE, info = cid)
    expect_true(all(g$noise_mode %in% c("measurement_homoscedastic", "measurement_heteroscedastic")),
                info = cid)
  }
})

test_that("measurement-noise reference: sample estimates by region", {
  skip_me()
  ref <- post(readRDS(ref_me_path())); res <- post(run_me())
  for (cid in EXPECTED_CIDS) {
    s <- res$plates[[cid]]$samples; r <- ref$plates[[cid]]$samples
    if (is.null(r)) { expect_null(s, info = cid); next }
    expect_equal(nrow(s), nrow(r), info = cid)
    in_range <- r$pcov_gate_class %in% "inDynamicRange"
    ok <- in_range & is.finite(s$final_concentration) & is.finite(r$final_concentration) &
          r$final_concentration > 0
    if (any(ok)) {
      rel <- abs(s$final_concentration[ok] - r$final_concentration[ok]) / r$final_concentration[ok]
      expect_true(all(rel < TOL_IN_RANGE_REL),
                  info = sprintf("curve %s in range: max relative diff %.4f", cid, max(rel)))
    }
    out <- !in_range & is.finite(s$predicted_log10_concentration) &
           is.finite(r$predicted_log10_concentration)
    if (any(out)) {
      d <- abs(s$predicted_log10_concentration[out] - r$predicted_log10_concentration[out])
      expect_true(all(d < TOL_OUT_LOG10),
                  info = sprintf("curve %s outside range: max log10 diff %.3f", cid, max(d)))
    }
    expect_true(any(in_range), info = paste("curve", cid, "has samples in the dynamic range"))
    same_gate <- mean(as.character(s$pcov_gate_class) == as.character(r$pcov_gate_class))
    expect_false(is.na(same_gate), info = cid)
    expect_gte(same_gate, MIN_CLASS_AGREE)
  }
})

test_that("grid predicted_concentration matches the June snapshot (its grid included noise)", {
  skip_me()
  skip_if_not(file.exists(ref_june_path()))
  june <- readRDS(ref_june_path()); res <- run_me()
  for (cid in EXPECTED_CIDS) {
    ref_pc <- june$plates[[cid]]$grid$predicted_concentration
    cur_pc <- res$plates[[cid]]$grid$predicted_concentration
    both_ok <- is.finite(ref_pc) & is.finite(cur_pc)
    if (any(both_ok))
      expect_equal(cur_pc[both_ok], ref_pc[both_ok], tolerance = 0.05,
                   info = paste("Predicted conc shifted for curve_id", cid))
    expect_equal(is.na(cur_pc), is.na(ref_pc), info = paste("NA pattern changed for curve_id", cid))
  }
})
