# =============================================================================
# data.R — Roxygen documentation for curveRbayes datasets
# =============================================================================


#' Simulated Bead-Based Immunoassay Example Dataset
#'
#' A synthetic multi-plate bead-based immunoassay dataset with two
#' antigens (`alpha`, `beta`) measured on three plates each, giving
#' six unique `curve_id` values.
#' Used in `vignette("bayesian-quickstart", package = "curveRbayes")`.
#'
#' @name bead_assay_example
#' @docType data
#' @keywords datasets
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{standards}{Data frame of calibration standards.
#'     Columns: `curve_id` (integer), `antigen` (character),
#'     `plate` (integer), `dilution` (numeric),
#'     `mfi` (numeric, raw median fluorescence intensity).}
#'   \item{samples}{Data frame of test samples.
#'     Columns: `curve_id`, `sampleid`, `dilution`, `mfi`.}
#'   \item{blanks}{Data frame of blank-well readings, or \code{NULL}
#'     if no blanks were recorded.
#'     Columns: `curve_id`, `mfi`.}
#'   \item{response_var}{Character. Name of the response column
#'     (`"mfi"`).}
#'   \item{indep_var}{Character. Name of the independent variable
#'     column used for concentration computation (`"dilution"`).}
#'   \item{curve_id_lookup}{Data frame mapping integer `curve_id`
#'     values to human-readable labels.
#'     Columns: `curve_id`, `antigen`, `plate`.}
#' }
#'
#' @details
#' Concentrations are simulated on a seven-point dilution series from
#' 0.001 to 30 AU/mL.  The undiluted standard concentration is 30
#' AU/mL (`std_curve_conc = 30`).  Both antigens share the same
#' dilution scheme but differ in signal level and curve shape to
#' exercise the inter-antigen variability path of the hierarchical
#' model.
#'
#' @source Simulated — no real patient data.
#'
#' @seealso \code{vignette("bayesian-quickstart", package = "curveRbayes")}
NULL


#' Simulated ELISA Example Dataset
#'
#' A synthetic single-plate ELISA dataset with one antigen measured
#' on a standard 96-well plate.  Intended as a minimal worked example
#' for users coming from a plate-reader workflow.
#'
#' @name elisa_assay_example
#' @docType data
#' @keywords datasets
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{standards}{Data frame of calibration standards.
#'     Columns: `curve_id` (integer), `dilution` (numeric),
#'     `od` (numeric, raw optical density).}
#'   \item{samples}{Data frame of test samples.
#'     Columns: `curve_id`, `sampleid`, `dilution`, `od`.}
#'   \item{blanks}{Data frame of blank-well readings, or \code{NULL}.
#'     Columns: `curve_id`, `od`.}
#'   \item{response_var}{Character. Name of the response column
#'     (`"od"`).}
#'   \item{indep_var}{Character. Name of the independent variable
#'     column (`"dilution"`).}
#'   \item{curve_id_lookup}{Data frame with columns `curve_id` and
#'     `label`.}
#' }
#'
#' @details
#' OD values are simulated on an eight-point dilution series.
#' The dataset is kept deliberately simple — one antigen, one plate —
#' so that the full preprocessing-to-back-calculation pipeline can be
#' demonstrated without the multi-plate complexity of
#' \code{\link{bead_assay_example}}.
#'
#' @source Simulated — no real patient data.
#'
#' @seealso \code{\link{bead_assay_example}}
NULL
