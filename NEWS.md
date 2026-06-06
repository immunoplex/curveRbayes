# curveRbayes 0.2.0

* Initial release.
* `fit_calibration_bayes()` — hierarchical Bayesian calibration across one
  or more `curve_id` values simultaneously via Stan. Returns a
  `calibration_result_multiplate` matching the curveRfreq output contract.
* Five Stan model families: `logistic4`, `logistic5`, `loglogistic4`,
  `loglogistic5`, `gompertz4` — all with non-centred parameterisation (NCP)
  and `reduce_sum` map-reduce likelihood for multi-core speedup.
* `compile_stan_model()` — compiles (and caches) a Stan model via cmdstanr.
* `fit_bayes_single()` — runs HMC/NUTS sampling for one model family across
  all curve_ids and returns posterior draws plus NUTS diagnostics
  (divergences, max treedepth, E-BFMI).
* `extract_curve_params()` — extracts per-curve posterior summaries
  (mean, SD, 2.5/50/97.5 quantiles) from a fitted model object.
* `build_stan_data()` — converts preprocessed stacked standards into the
  named-list format expected by the Stan `data {}` block, including all
  data-adaptive prior scalars.
* `compute_dynamic_priors()` — constructs weakly informative, data-adaptive
  hyperpriors from the preprocessed standards. Supports `fixed_a` soft
  constraint for poorly-identified lower asymptotes.
* `predict_grid_bayes()` — CDAN (Concentration-Dependent Assay Noise)
  precision grid: three-step procedure of posterior draw → forward
  evaluation → Student-t noise injection → analytical back-calculation,
  producing `pcov` and `pcov_rmse` precision profiles at every grid point.
* `predict_samples_bayes()` — back-calculates test-sample concentrations
  from posterior draws without noise injection (observed response is already
  the noisy measurement). Returns `final_concentration`, `se_concentration`,
  `pcov`, and `pcov_pass` per sample.
* `compute_loo()` / `compare_models_loo()` — PSIS-LOO cross-validation and
  Bayesian stacking weights via the `loo` package.
* `summary_table_bayes()` / `collect_samples_bayes()` — tidy extraction from
  `calibration_result_multiplate` objects (multiplate and legacy single-curve
  formats both supported).
* Two eligibility gates active on the Bayesian path: `rel_se` (posterior
  SD / |mean| per parameter) and `dynamic_range` (log10 upper/lower
  asymptote ratio). `at_bound` and `vcov_condition` gates are bypassed
  (no hard constraints; no vcov matrix).
* Global eligibility: a model must pass both gates on **all** curve_ids
  to be eligible for back-calculation.
* `bead_assay_example` synthetic dataset: two antigens × three plates,
  six `curve_id` values.
