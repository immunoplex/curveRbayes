# Package index

## Main Entry Point

The function most users call directly. Expects preprocessed standards on
the fitting scale (use
[`curveRcore::preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.html)
upstream). Fits all `curve_id` values simultaneously via hierarchical
Stan models and returns a `calibration_result_multiplate`.

- [`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
  : Fit Bayesian Hierarchical Calibration Curves

## Multi-Curve Helpers

Convenience extractors that operate on `calibration_result_multiplate`
objects returned by
[`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md).
Both functions also accept the legacy single-`calibration_result`
format.

- [`summary_table_bayes()`](https://immunoplex.github.io/curveRbayes/reference/summary_table_bayes.md)
  : Extract a Per-Curve Summary Table from Bayesian Results
- [`collect_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/collect_samples_bayes.md)
  : Collect All Sample Predictions from Bayesian Results
- [`collect_standards_bayes()`](https://immunoplex.github.io/curveRbayes/reference/collect_standards_bayes.md)
  : Collect All Standard Data from Bayesian Results
- [`collect_blanks_bayes()`](https://immunoplex.github.io/curveRbayes/reference/collect_blanks_bayes.md)
  : Collect All Blank Data from Bayesian Results

## Stan Compilation and MCMC Fitting

Lower-level functions that compile Stan models and run HMC/NUTS
sampling. Called internally by
[`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
but exported for users who need fine-grained control over model
compilation, sampling parameters, or multi-step workflows.

- [`compile_stan_model()`](https://immunoplex.github.io/curveRbayes/reference/compile_stan_model.md)
  : Compile a curveRbayes Stan Model (Cached)
- [`fit_bayes_single()`](https://immunoplex.github.io/curveRbayes/reference/fit_bayes_single.md)
  : Fit a Single Model Family via MCMC
- [`extract_curve_params()`](https://immunoplex.github.io/curveRbayes/reference/extract_curve_params.md)
  : Extract Curve-Level Posterior Summaries

## Stan Data and Priors

Functions that prepare inputs for Stan.
[`compute_dynamic_priors()`](https://immunoplex.github.io/curveRbayes/reference/compute_dynamic_priors.md)
derives weakly informative hyperpriors from the preprocessed data range.
[`build_stan_data()`](https://immunoplex.github.io/curveRbayes/reference/build_stan_data.md)
assembles all Stan inputs — observations, curve indices, and prior
scalars — into the named list expected by the Stan `data {}` block.

- [`compute_dynamic_priors()`](https://immunoplex.github.io/curveRbayes/reference/compute_dynamic_priors.md)
  : Compute Data-Adaptive Priors for Stan Models
- [`build_stan_data()`](https://immunoplex.github.io/curveRbayes/reference/build_stan_data.md)
  : Build Stan Data List for a Model Family

## LOO-CV Model Selection

PSIS-LOO cross-validation and Bayesian stacking weights.
[`compute_loo()`](https://immunoplex.github.io/curveRbayes/reference/compute_loo.md)
extracts the `log_lik` generated quantity and computes a `loo` object.
[`compare_models_loo()`](https://immunoplex.github.io/curveRbayes/reference/compare_models_loo.md)
runs LOO for every fitted model and returns the `loo_compare()` table
plus stacking weights. Called automatically by
[`fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.md)
when more than one model is fitted.

- [`compute_loo()`](https://immunoplex.github.io/curveRbayes/reference/compute_loo.md)
  : Compute LOO-CV for a Fitted Bayesian Model
- [`compare_models_loo()`](https://immunoplex.github.io/curveRbayes/reference/compare_models_loo.md)
  : Compare Models via LOO-CV and Stacking Weights

## Posterior Prediction and CDAN Precision

Posterior predictive grid construction and test-sample back-calculation.
[`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
implements the three-step CDAN procedure (posterior draw, forward
evaluation, Student-t noise injection, analytical inversion) to produce
a full precision profile.
[`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
back-calculates observed test-sample responses without noise injection,
since the observed response is already the noisy measurement.

- [`predict_grid_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_grid_bayes.md)
  : Predict Grid Response from Posterior Draws (Bayesian)
- [`predict_samples_bayes()`](https://immunoplex.github.io/curveRbayes/reference/predict_samples_bayes.md)
  : Back-Calculate Sample Concentrations from Posterior Draws

## Example Datasets

Synthetic immunoassay datasets for testing and documentation.

- [`bead_assay_example`](https://immunoplex.github.io/curveRbayes/reference/bead_assay_example.md)
  : Simulated Bead-Based Immunoassay Example Dataset
- [`elisa_assay_example`](https://immunoplex.github.io/curveRbayes/reference/elisa_assay_example.md)
  : Simulated ELISA Example Dataset
