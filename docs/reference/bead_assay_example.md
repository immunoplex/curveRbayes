# Simulated Bead-Based Immunoassay Example Dataset

A synthetic multi-plate bead-based immunoassay dataset with two antigens
(`alpha`, `beta`) measured on three plates each, giving six unique
`curve_id` values. Used in
[`vignette("bayesian-quickstart", package = "curveRbayes")`](https://immunoplex.github.io/curveRbayes/articles/bayesian-quickstart.md).

## Format

A named list with the following elements:

- standards:

  Data frame of calibration standards. Columns: `curve_id` (integer),
  `antigen` (character), `plate` (integer), `dilution` (numeric), `mfi`
  (numeric, raw median fluorescence intensity).

- samples:

  Data frame of test samples. Columns: `curve_id`, `sampleid`,
  `dilution`, `mfi`.

- blanks:

  Data frame of blank-well readings, or `NULL` if no blanks were
  recorded. Columns: `curve_id`, `mfi`.

- response_var:

  Character. Name of the response column (`"mfi"`).

- indep_var:

  Character. Name of the independent variable column used for
  concentration computation (`"dilution"`).

- curve_id_lookup:

  Data frame mapping integer `curve_id` values to human-readable labels.
  Columns: `curve_id`, `antigen`, `plate`.

## Source

Simulated — no real patient data.

## Details

Concentrations are simulated on a seven-point dilution series from 0.001
to 30 AU/mL. The undiluted standard concentration is 30 AU/mL
(`std_curve_conc = 30`). Both antigens share the same dilution scheme
but differ in signal level and curve shape to exercise the inter-antigen
variability path of the hierarchical model.

## See also

[`vignette("bayesian-quickstart", package = "curveRbayes")`](https://immunoplex.github.io/curveRbayes/articles/bayesian-quickstart.md)
