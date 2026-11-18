#!/usr/bin/env Rscript
# =============================================================================
# generate_reference_bayes_me.R
#
# Frozen reference for test-reference-bayes-me.R: the same fit as
# generate_reference_bayes.R (curve_ids 1-3 of bead_assay_example, same seed),
# but with the DEFAULT back-calculation: measurement noise included
# (include_measurement_error = TRUE) and seeded (curveRbayes >= 0.5.2), so the
# noise draws are reproducible.
#
# Run once from the curveRbayes root, and again only after an intentional
# change to the fit or the back-calculation (say why in NEWS.md):
#
#   Rscript data-raw/generate_reference_bayes_me.R
#
# Saved to tests/testthat/fixtures/reference_bayes_curve123_me.rds, with the
# versions it was made under in attr(, "generated_with").
# =============================================================================
library(curveRcore)
devtools::load_all(quiet = TRUE)   # the curveRbayes source being released

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
  run_loo = TRUE, n_draws_predict = 500,
  include_measurement_error = TRUE,
  verbose = FALSE
)
attr(ref, "generated_with") <- list(
  date = format(Sys.Date()),
  curveRbayes = as.character(utils::packageVersion("curveRbayes")),
  curveRcore = as.character(utils::packageVersion("curveRcore")),
  cmdstan = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA),
  r = paste(R.version$major, R.version$minor, sep = "."))

saveRDS(ref, "tests/testthat/fixtures/reference_bayes_curve123_me.rds")
cat("Saved reference_bayes_curve123_me.rds\n")
str(attr(ref, "generated_with"))
