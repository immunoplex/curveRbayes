# Censoring, predictive intervals and threshold probabilities for test
# samples, including responses beyond the curve's asymptotes (previously
# untested). Minimal sampler settings: correctness of the bookkeeping, not
# of the posterior.

library(curveRcore)

skip_if_no_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) skip("cmdstanr not available")
  tryCatch(cmdstanr::cmdstan_path(), error = function(e) skip("CmdStan not installed"))
}

data("bead_assay_example", package = "curveRcore")
oor_std <- curveRcore::preprocess_standards(
  data = bead_assay_example$standards[bead_assay_example$standards$curve_id == 1, ],
  antigen_settings = list(standard_curve_concentration = 10000),
  response_variable = "mfi", independent_variable = "concentration",
  is_log_response = TRUE, is_log_independent = TRUE, apply_prozone = TRUE)$data
oor_samp <- bead_assay_example$samples[bead_assay_example$samples$curve_id == 1, ]
oor_samp <- oor_samp[seq_len(min(6, nrow(oor_samp))), ]
extremes <- oor_samp[1:2, ]
extremes$sampleid <- c("very_low", "very_high")
extremes$mfi <- c(1e-3, 1e9)
oor_samp <- rbind(oor_samp, extremes)

new_cols <- c("conc_q_lo", "conc_q_med", "conc_q_hi", "p_below_lloq", "p_above_uloq",
              "p_below_cutoff", "frac_draws_below_a", "frac_draws_above_d")

test_that("Bayesian samples carry intervals, probabilities and censoring", {
  skip_on_cran()
  skip_if_no_cmdstan()
  mp <- fit_calibration_bayes(
    standards = oor_std, samples = oor_samp, response_var = "mfi",
    model_names = "logistic4", std_curve_conc = 10000,
    chains = 1L, warmup = 200L, sampling = 200L, seed = 1L,
    n_draws_predict = 200L, decision_cutoff = 50, run_loo = FALSE)
  mp <- curveRcore::compute_detection_limits_multiplate(mp)
  mp <- curveRcore::classify_pcov_gates_multiplate(mp)
  mp <- curveRcore::classify_censoring_multiplate(mp)
  s  <- mp$plates[[1]]$samples

  expect_true(all(c(new_cols, "censor_class", "censor_type") %in% names(s)))
  expect_equal(mp$meta$decision_cutoff, 50)

  lo <- s[s$sampleid == "very_low", ]; hi <- s[s$sampleid == "very_high", ]
  expect_equal(lo$frac_draws_below_a, 1)
  expect_equal(lo$conc_q_hi, 0)
  expect_equal(lo$p_below_cutoff, 1)
  expect_equal(lo$censor_type, "left")
  expect_equal(hi$frac_draws_above_d, 1)
  expect_true(is.na(hi$conc_q_lo))
  expect_equal(hi$censor_type, "right")

  ok <- is.finite(s$conc_q_lo) & is.finite(s$conc_q_hi)
  expect_true(all(s$conc_q_lo[ok] <= s$conc_q_med[ok] & s$conc_q_med[ok] <= s$conc_q_hi[ok]))
})
