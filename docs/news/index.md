# Changelog

## curveRbayes 0.5.1

- Maintenance release: no functional changes. Rebuilt and documented
  against curveRcore 0.5.2 (whose `getting-started` vignette now covers
  censoring and precision weights).

## curveRbayes 0.5.0

- **Predictive intervals and threshold probabilities for test samples.**
  [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
  now also inverts every posterior draw with out-of-range responses kept
  as censored
  ([`curveRcore::invert_with_bounds()`](https://immunoplex.github.io/curveRcore/reference/invert_with_bounds.html))
  and adds `conc_q_lo` / `conc_q_med` / `conc_q_hi`, `p_below_lloq`,
  `p_above_uloq`, `p_below_cutoff`, `frac_draws_below_a` and
  `frac_draws_above_d` via
  [`curveRcore::summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.html).
  New `decision_cutoff` argument on
  [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  (recorded in `$meta`).
- `se_concentration` / `pcov` / `final_concentration` are byte-identical
  to 0.4.3 for the same posterior and seed: the new path reuses each
  draw’s already-sampled response and adds no RNG calls. The new
  `frac_draws_*` columns show how many draws those legacy summaries
  drop.
- Requires `curveRcore (>= 0.5.0)`.
- Vignette: replaced the stale “noise is NOT injected for test samples”
  section — samples follow `include_measurement_error` (default `TRUE`),
  the same setting as the grid.
- **Bug fix:**
  [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  now passes `is_log_response` to
  [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
  (the argument had been commented out), so test samples are
  back-calculated on the same response scale as the standards.
  Previously samples always assumed a log10 response; results change
  only for fits with `is_log_response = FALSE` (set through the settings
  cascade).
- Known, pre-existing: `test-reference-bayes.R` “sample predictions
  match reference” fails on the stored fixture with the 0.4.3 code as
  well; regenerate the fixture.

## curveRbayes 0.4.3 (2026-09-26)

- Declared the real dependency floor: `Imports: curveRcore (>= 0.4.0)`.
  [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  passes `population =` to
  [`new_calibration_result_multiplate()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result_multiplate.html),
  which curveRcore only accepts from 0.4.0. With no floor, a clean
  install could resolve curveRcore 0.2.0 and then fail every fit with
  `unused argument (population = population)` at the end of the run.
- `hierarchical_loglogistic4.stan` now uses
  `sigma_log_c ~ cauchy(0, 0.5)`, completing the half-Cauchy switch from
  0.4.2+. That family estimates the location on the log scale and so
  names its between-plate SD `sigma_log_c`; the original patch matched
  `sigma_(a|d|c|log_b)` and skipped it, leaving one family with a
  half-Normal location tau while the other four were half-Cauchy. The
  scale is 0.5, not the previous 1.0, because a half-Cauchy(0, 1) on a
  log-EC50 SD admits implausibly large between-plate spread.
- Version bump so the tag uniquely identifies the code. `main` had
  carried the LOO-stacking fix, the Gompertz overflow guard, the
  `reduce_sum` threading and the half-Cauchy taus while still declaring
  0.4.2, so two materially different package states both answered
  `packageVersion("curveRbayes")` with “0.4.2”.

## curveRbayes 0.4.2

- Created a new pcov_gate_class and changed the basis for pcov_pass
  classifcations.

## curveRbayes 0.4.1

- fit_calibration_bayes() now passes pcov_threshold explicitly into
  predict_samples_bayes() instead of relying on that function’s own
  default, so per-sample pcov_pass always reflects the study’s
  configured threshold.

## curveRbayes 0.4.0 (2026-07-29)

- **Unified the precision-profile variance definition** across the grid
  and the per-sample back-calculation via a new
  `include_measurement_error` argument to
  [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  (default `TRUE`), threaded identically into
  [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
  and
  [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md).
  This closes the gap where Bayesian sample points floated ~9.6 %CV
  below the precision profile with a concentration-dependent slope: with
  the default they now lie on the measurement/CDAN profile.
- [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
  and
  [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
  now share a single
  [`.obs_noise_sigma()`](https://immunoplex.github.io/curveRbayes/reference/dot-obs_noise_sigma.md)
  noise definition, so the grid and sample paths cannot diverge again.
- `include_measurement_error = FALSE` yields a
  curve/parameter-uncertainty-only profile; in that mode the grid
  inverts a fixed posterior-mean reference response across draws
  (non-degenerate) to match the sample path.
- The active mode is recorded in the `noise_mode` column and in the
  result `meta`.

## curveRbayes 0.2.0

- Initial release.
- [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  — hierarchical Bayesian calibration across one or more `curve_id`
  values simultaneously via Stan. Returns a
  `calibration_result_multiplate` matching the curveRfreq output
  contract.
- Five Stan model families: `logistic4`, `logistic5`, `loglogistic4`,
  `loglogistic5`, `gompertz4` — all with non-centred parameterisation
  (NCP) and `reduce_sum` map-reduce likelihood for multi-core speedup.
- [`compile_stan_model()`](https://immunoplex.github.io/curveRbayes/reference/compile_stan_model.md)
  — compiles (and caches) a Stan model via cmdstanr.
- [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md)
  — runs HMC/NUTS sampling for one model family across all curve_ids and
  returns posterior draws plus NUTS diagnostics (divergences, max
  treedepth, E-BFMI).
- [`extract_curve_params()`](https://immunoplex.github.io/curveRbayes/reference/extract_curve_params.md)
  — extracts per-curve posterior summaries (mean, SD, 2.5/50/97.5
  quantiles) from a fitted model object.
- [`build_stan_data()`](https://immunoplex.github.io/curveRbayes/reference/build_stan_data.md)
  — converts preprocessed stacked standards into the named-list format
  expected by the Stan `data {}` block, including all data-adaptive
  prior scalars.
- [`compute_dynamic_priors()`](https://immunoplex.github.io/curveRbayes/reference/compute_dynamic_priors.md)
  — constructs weakly informative, data-adaptive hyperpriors from the
  preprocessed standards. Supports `fixed_a` soft constraint for
  poorly-identified lower asymptotes.
- [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
  — CDAN (Concentration-Dependent Assay Noise) precision grid:
  three-step procedure of posterior draw → forward evaluation →
  Student-t noise injection → analytical back-calculation, producing
  `pcov` and `pcov_rmse` precision profiles at every grid point.
- [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
  — back-calculates test-sample concentrations from posterior draws
  without noise injection (observed response is already the noisy
  measurement). Returns `final_concentration`, `se_concentration`,
  `pcov`, and `pcov_pass` per sample.
- [`compute_loo()`](https://immunoplex.github.io/curveRbayes/reference/compute_loo.md)
  /
  [`compare_models_loo()`](https://immunoplex.github.io/curveRbayes/reference/compare_models_loo.md)
  — PSIS-LOO cross-validation and Bayesian stacking weights via the
  `loo` package.
- [`summary_table_bayes()`](https://immunoplex.github.io/curveRbayes/reference/summary_table_bayes.md)
  /
  [`collect_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/collect_samples_bayes.md)
  — tidy extraction from `calibration_result_multiplate` objects
  (multiplate and legacy single-curve formats both supported).
- Two eligibility gates active on the Bayesian path: `rel_se` (posterior
  SD / \|mean\| per parameter) and `dynamic_range` (log10 upper/lower
  asymptote ratio). `at_bound` and `vcov_condition` gates are bypassed
  (no hard constraints; no vcov matrix).
- Global eligibility: a model must pass both gates on **all** curve_ids
  to be eligible for back-calculation.
- `bead_assay_example` synthetic dataset: two antigens × three plates,
  six `curve_id` values.
