#!/usr/bin/env Rscript
# =============================================================================
# generate_reference_bayes.R
#
# Generates the frozen reference RDS for the Bayesian regression test.
# Run this ONCE from the curveRbayes root with:
#
#   Rscript data-raw/generate_reference_bayes.R
#
# The reference is saved to:
#   tests/testthat/fixtures/reference_bayes_bead_alpha.rds
#
# Requirements: curveRcore, curveRbayes, cmdstanr, CmdStan installed.
# =============================================================================
library(curveRcore)
library(curveRbayes)

data("bead_assay_example", package = "curveRcore")
std_raw <- bead_assay_example$standards[
  bead_assay_example$standards$curve_id %in% c(1, 2, 3), ]
prepped <- preprocess_standards(
  data = std_raw,
  antigen_settings = list(standard_curve_concentration = 10000),
  response_variable = "mfi", independent_variable = "concentration",
  is_log_response = TRUE, is_log_independent = TRUE, apply_prozone = TRUE
)
samples <- bead_assay_example$samples[
  bead_assay_example$samples$curve_id %in% c(1, 2, 3), ]

ref <- fit_calibration_bayes(
  standards = prepped$data, samples = samples, response_var = "mfi",
  model_names = c("logistic4", "gompertz4"),
  is_log_response = TRUE, is_log_independent = TRUE, std_curve_conc = 10000,
  chains = 4, warmup = 1000, sampling = 1000,
  adapt_delta = 0.9, seed = 20260529,
  run_loo = TRUE, n_draws_predict = 500
)

saveRDS(ref, "tests/testthat/fixtures/reference_bayes_curve123.rds")
cat("Saved reference_bayes_curve123.rds\n")
