# curveRbayes

Bayesian hierarchical calibration curves for the curveR suite.

## Installation

```r
# Install curveRcore first
devtools::install_github("immunoplex/curveRcore")
# Then curveRbayes
devtools::install_github("immunoplex/curveRbayes")
```

### Stan toolchain

curveRbayes fits models via **cmdstanr**.
Install it once after installing the package:

```r
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           getOption("repos")))
cmdstanr::install_cmdstan()   # downloads and builds CmdStan
```

## Quick start

```r
library(curveRbayes)

data(bead_assay_example, package = "curveRbayes")

# Preprocess standards (upstream of fitting)
std_pre <- curveRcore::preprocess_standards(
  data                 = bead_assay_example$standards,
  antigen_settings     = antigen_settings,
  response_variable    = bead_assay_example$response_var,
  independent_variable = bead_assay_example$indep_var,
  is_log_response      = TRUE,
  is_log_independent   = TRUE
)

# Fit all curves simultaneously
mp <- fit_calibration_bayes(
  standards          = std_pre$data,
  samples            = bead_assay_example$samples,
  response_var       = "mfi",
  model_names        = c("logistic4", "gompertz4"),
  std_curve_conc     = 30,
  chains             = 4L,
  warmup             = 1000L,
  sampling           = 1000L,
  seed               = 42
)

# One-row-per-curve summary
summary_table_bayes(mp)

# All sample back-calculations
collect_samples_bayes(mp)
```

---

## Quick Reference

| Task | Command |
|---|---|
| Regenerate man pages | `devtools::document()` |
| Full site rebuild | `pkgdown::build_site(dest_dir = "docs")` |
| Reference pages only | `pkgdown::build_reference()` |
| Vignette only | `pkgdown::build_articles()` |
| Home page only | `pkgdown::build_home()` |
| Check no topics are unassigned | `setdiff(pkg$topics$name, unlist(lapply(cfg$reference, "[[", "contents")))` |
| Preview locally | Open `docs/index.html` |
| Push and deploy | `git add docs/ && git commit -m "..." && git push` |
| Check CmdStan installation | `cmdstanr::cmdstan_version()` |
| Recompile a Stan model | `compile_stan_model("logistic4")` |

---

## Stan models

curveRbayes ships five pre-written Stan models in `inst/stan/`.
All use a **non-centred parameterisation** (NCP) to eliminate
Neal's funnel geometry and a **`reduce_sum`** map-reduce likelihood
for multi-core speedup.

| File | Model family | Parameters |
|---|---|---|
| `hierarchical_logistic4.stan` | 4-parameter logistic | A, B, C, D |
| `hierarchical_logistic5.stan` | 5-parameter logistic | A, B, C, D, G |
| `hierarchical_loglogistic4.stan` | 4-parameter log-logistic | A, B, C, D |
| `hierarchical_loglogistic5.stan` | 5-parameter log-logistic | A, B, C, D, G |
| `hierarchical_gompertz4.stan` | 4-parameter Gompertz | A, B, C, D |

When `is_log_independent = TRUE` and `is_log_response = TRUE`, `loglogistic4`
is mathematically equivalent to `logistic4` and is automatically dropped
from the candidate set.

See `vignette("bayesian-quickstart", package = "curveRbayes")` for the
full worked example including Stan model internals, MCMC diagnostics,
LOO-CV model selection, and the CDAN precision profiling procedure.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `curveRcore` not found during CI | Add `install_github("immunoplex/curveRcore")` before `setup-r-dependencies` in the workflow |
| `cmdstanr` not found | Install from Stan universe: `install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))` |
| CmdStan not installed | Run `cmdstanr::install_cmdstan()` |
| Stan file not found (`mustWork = TRUE` error) | Run `devtools::install()` so `inst/stan/` is copied to the package library |
| Divergent transitions > 10 | Increase `adapt_delta` toward 0.99; inspect pairs plots with `bayesplot::mcmc_pairs()` |
| Low E-BFMI (< 0.2) | Check prior–likelihood conflict; consider supplying `fixed_a` for poorly identified lower asymptote |
| Max treedepth hits frequent | Increase `max_treedepth` in `fit_bayes_single()` or widen priors |
| Vignette fails during `R CMD check` | Add `eval = requireNamespace("curveRcore", quietly = TRUE)` to the setup chunk |
| `bead_assay_example` not found | Run `usethis::use_data(bead_assay_example)` and add a `R/data.R` roxygen block |
| Function missing from Reference | Check `@export` is present; re-run `devtools::document()` |
| LOO Pareto-k > 0.7 | That calibrator is high-leverage; consider moment-matching via `loo::loo(moment_match = TRUE)` |
| `docs/` not served by GitHub | Settings → Pages: branch `main`, folder `/docs` |

---

See `vignette("bayesian-quickstart", package = "curveRbayes")` for
the full worked example.
