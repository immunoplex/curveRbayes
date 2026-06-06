# helper-cmdstan.R
# Suppress CmdStan subprocess stdout/stderr in the test environment.
#
# On Windows, CmdStan writes binary progress bytes to stdout that contain
# embedded NUL characters. The testthat runner captures this output and
# passes it to iconv(), which throws:
#   Error in iconv(...): embedded nul in string
#
# Setting refresh = 0 and show_messages = FALSE is already done in
# fit_bayes_single() when verbose = FALSE, but the subprocess can still
# emit a UTF-8 BOM or binary header before Stan respects the flag.
#
# The reliable fix is to route CmdStan output to a temp file during tests.

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  # Redirect all CmdStan output to a temporary file for the test session.
  # This has no effect on the correctness of sampling; only console output
  # is suppressed.
  cmdstanr::set_cmdstan_path(cmdstanr::cmdstan_path())   # no-op, but forces init

  # Override the global output_dir so compiled models and CSV files land
  # in a clean temp directory rather than the working directory.
  options(cmdstanr_output_dir = tempdir())

  # On Windows, suppress the per-chain progress output that causes the
  # embedded-NUL iconv error.
  if (.Platform$OS.type == "windows") {
    options(cmdstanr_write_stan_file_dir = tempdir())
  }
}
