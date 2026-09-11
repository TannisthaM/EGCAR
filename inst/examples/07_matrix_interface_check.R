# Run after a clean source reinstall and restart of R.
# RStudio: source(system.file("examples", "07_matrix_interface_check.R",
#                            package = "egcar", mustWork = TRUE))
# Defaults: signal 0.8, ranks 1/2/5, five folds, five CV workers, small grids.
# A serial-only diagnostic: Sys.setenv(EGCAR_MATRIX_CHECK_WORKERS = "1")
if (!requireNamespace("egcar", quietly = TRUE)) stop("Install egcar first.")
if (utils::packageVersion("egcar") < "0.2.4")
  stop("This check requires the patched egcar 0.2.2 package. Reinstall and restart R.")
local({
  workers <- suppressWarnings(as.integer(Sys.getenv("EGCAR_MATRIX_CHECK_WORKERS", "5")))
  if (is.na(workers) || workers < 1L || workers > 5L)
    stop("EGCAR_MATRIX_CHECK_WORKERS must be an integer from 1 to 5.")
  checks <- new.env(parent = globalenv())
  sys.source(system.file("validation", "matrix_interface.R", package = "egcar",
                          mustWork = TRUE), envir = checks)
  report <- checks$run_egcar_matrix_interface_check(workers = workers)
  print(report$native, row.names = FALSE)
  print(report$fits, row.names = FALSE)
  print(report$cv, row.names = FALSE)
  out <- Sys.getenv("EGCAR_MATRIX_CHECK_OUTPUT", "egcar_matrix_interface_check")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(report, file.path(out, "matrix_interface_report.rds"))
  utils::write.csv(report$native, file.path(out, "native_matrix_checks.csv"), row.names = FALSE)
  utils::write.csv(report$fits, file.path(out, "fit_matrix_checks.csv"), row.names = FALSE)
  utils::write.csv(report$cv, file.path(out, "cv_matrix_checks.csv"), row.names = FALSE)
  utils::capture.output(report$session_info, file = file.path(out, "sessionInfo.txt"))
  message("Check reports saved in ", normalizePath(out))
})
