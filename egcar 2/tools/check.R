#!/usr/bin/env Rscript
# Run in the package root: Rscript tools/check.R .
# Optional comparisons: EGCAR_TEST_OPTIONAL=true Rscript tools/check.R .
.run_checks <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  path <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
  description <- file.path(path, "DESCRIPTION")
  if (!file.exists(description)) stop("Not an R package root: ", path)
  required <- c("Rcpp", "RcppArmadillo", "testthat")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the check dependencies first: ", paste(missing, collapse = ", "))
  for (dir in c("R", "inst/examples", "inst/validation", "tests")) {
    for (f in list.files(file.path(path, dir), pattern = "\\.[Rr]$",
                         full.names = TRUE, recursive = TRUE)) parse(f)
  }
  for (f in list.files(file.path(path, "man"), pattern = "\\.Rd$", full.names = TRUE))
    tools::parse_Rd(f)
  desc <- read.dcf(description)
  archive <- paste0(desc[1L, "Package"], "_", desc[1L, "Version"], ".tar.gz")
  R <- file.path(R.home("bin"), "R")
  oldwd <- setwd(dirname(path)); on.exit(setwd(oldwd), add = TRUE)
  old_suggests <- Sys.getenv("_R_CHECK_FORCE_SUGGESTS_", unset = NA_character_)
  on.exit({
    if (is.na(old_suggests)) Sys.unsetenv("_R_CHECK_FORCE_SUGGESTS_") else
      Sys.setenv(`_R_CHECK_FORCE_SUGGESTS_` = old_suggests)
  }, add = TRUE)
  status <- system2(R, c("CMD", "build", "--no-build-vignettes", shQuote(path)))
  if (status != 0L) stop("R CMD build failed.")
  if (!file.exists(archive)) stop("Build did not create expected archive: ", archive)
  Sys.setenv(`_R_CHECK_FORCE_SUGGESTS_` = "false")
  status <- system2(R, c("CMD", "check", "--no-manual", shQuote(archive)))
  if (status != 0L) stop("R CMD check reported a failure; inspect ",
                         desc[1L, "Package"], ".Rcheck/00check.log.")
}
.run_checks()
