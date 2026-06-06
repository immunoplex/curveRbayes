# Simulated ELISA Example Dataset

A synthetic single-plate ELISA dataset with one antigen measured on a
standard 96-well plate. Intended as a minimal worked example for users
coming from a plate-reader workflow.

## Format

A named list with the following elements:

- standards:

  Data frame of calibration standards. Columns: `curve_id` (integer),
  `dilution` (numeric), `od` (numeric, raw optical density).

- samples:

  Data frame of test samples. Columns: `curve_id`, `sampleid`,
  `dilution`, `od`.

- blanks:

  Data frame of blank-well readings, or `NULL`. Columns: `curve_id`,
  `od`.

- response_var:

  Character. Name of the response column (`"od"`).

- indep_var:

  Character. Name of the independent variable column (`"dilution"`).

- curve_id_lookup:

  Data frame with columns `curve_id` and `label`.

## Source

Simulated — no real patient data.

## Details

OD values are simulated on an eight-point dilution series. The dataset
is kept deliberately simple — one antigen, one plate — so that the full
preprocessing-to-back-calculation pipeline can be demonstrated without
the multi-plate complexity of
[`bead_assay_example`](https://immunoplex.github.io/curveRbayes/reference/bead_assay_example.md).

## See also

[`bead_assay_example`](https://immunoplex.github.io/curveRbayes/reference/bead_assay_example.md)
