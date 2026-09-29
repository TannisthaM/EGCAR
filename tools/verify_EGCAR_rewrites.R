#!/usr/bin/env Rscript
# Rscript tools/verify_EGCAR_rewrites.R [package_source_directory]
# The duplicate standalone solver was retired in 0.2.16. Validate the maintained
# implementation against independent mathematical references instead.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1L]] else "."
if (!requireNamespace("egcar",quietly=TRUE) || utils::packageVersion("egcar") < "0.2.16")
  stop("Install egcar >= 0.2.16 first.")
if (!requireNamespace("testthat",quietly=TRUE)) stop("Install testthat first.")
loadNamespace("egcar")
testthat::test_dir(file.path(root,"tests","testthat"),
  filter="sgca-paper-0216|efficiency-timing-0215|revision-020",
  reporter="summary",env=new.env(parent=asNamespace("egcar")),stop_on_failure=TRUE)
