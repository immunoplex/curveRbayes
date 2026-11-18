# Seeded measurement-noise draws (0.5.2): reproducible sample and grid
# estimates, independent of the other samples, caller's RNG untouched.
# Uses a synthetic posterior, so no CmdStan is needed.

fake_fit <- function(S = 400, hetero = FALSE) {
  set.seed(1)
  d <- list(`a[1]` = rnorm(S, 1.6, 0.05), `b[1]` = rnorm(S, 0.5, 0.02),
            `c_par[1]` = rnorm(S, 1.9, 0.04), `d[1]` = rnorm(S, 4.3, 0.03),
            sigma_obs = abs(rnorm(S, 0.05, 0.005)), nu = runif(S, 5, 30))
  if (hetero) { d$log_sigma0 <- rnorm(S, -3, 0.1); d$log_sigma_slope <- rnorm(S, 0.2, 0.02) }
  list(model_family = "logistic4", draws = d,
       stan_data = list(use_heteroscedastic_noise = as.integer(hetero)))
}
smp <- data.frame(sampleid = paste0("s", 1:6), mfi = 10^c(1.8, 2.4, 3.0, 3.6, 4.1, 4.25))
pred <- function(samples, seed, n_draws = 300, fit = fake_fit())
  predict_samples_bayes(samples, fit, response_variable = "mfi", n_draws = n_draws,
                        include_measurement_error = TRUE, seed = seed)
cols <- c("predicted_log10_concentration", "se_concentration", "pcov",
          "conc_q_lo", "conc_q_med", "conc_q_hi", "frac_draws_above_d")

test_that("the same seed gives identical sample estimates", {
  a <- pred(smp, seed = 42); b <- pred(smp, seed = 42)
  expect_identical(a[cols], b[cols])
  expect_false(isTRUE(all.equal(a$se_concentration, pred(smp, seed = 43)$se_concentration)))
})

test_that("the caller's random-number state is left untouched", {
  f <- fake_fit()                      # (builds the synthetic posterior with set.seed)
  set.seed(7); before <- .Random.seed
  invisible(pred(smp, seed = 42, fit = f))
  expect_identical(.Random.seed, before)
})

test_that("a sample's estimate does not depend on the other samples or their order", {
  full <- pred(smp, seed = 42)
  one  <- pred(smp[3, , drop = FALSE], seed = 42)
  expect_equal(one[cols], full[3, cols], ignore_attr = TRUE)
  rev  <- pred(smp[6:1, ], seed = 42)
  expect_equal(rev[6:1, cols], full[cols], ignore_attr = TRUE)
})

test_that("the grid is reproducible with a seed, and works without one", {
  g <- curveRcore::generate_prediction_grid(std_curve_conc = 10000, n_grid = 25,
                                            grid_min_conc = 1e-4, grid_max_conc = NULL,
                                            is_log_independent = TRUE)
  g1 <- predict_grid_bayes(g, fake_fit(), n_draws = 300, seed = 9)
  g2 <- predict_grid_bayes(g, fake_fit(), n_draws = 300, seed = 9)
  expect_identical(g1$se_concentration, g2$se_concentration)
  g0 <- predict_grid_bayes(g, fake_fit(), n_draws = 300)
  expect_true(all(is.finite(g0$pcov)))
  expect_true(all(is.finite(pred(smp, seed = NULL)$pcov)))
})

test_that("heteroscedastic noise is seeded too", {
  f <- fake_fit(hetero = TRUE)
  expect_identical(pred(smp, 5, fit = f)[cols], pred(smp, 5, fit = f)[cols])
})

test_that("the noise keeps its t(nu) distribution", {
  ns <- .noise_stream(S = 20000, n_draws = NULL, nu_draws = rep(8, 20000), seed = 3, curve_idx = 1)
  expect_equal(stats::var(ns$e), 8 / 6, tolerance = 0.05)    # var of t(8) = 8/6
  ns2 <- .noise_stream(S = 20000, n_draws = NULL, nu_draws = rep(8, 20000), seed = 3, curve_idx = 2)
  expect_false(identical(ns$e, ns2$e))                       # curves get their own stream
})
