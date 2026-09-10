#!/usr/bin/env Rscript
# Explicit installation helper; package fitting never installs anything.
args <- commandArgs(trailingOnly = TRUE)
packages <- c("Rcpp", "RcppArmadillo", "future", "future.apply", "RSpectra", "RhpcBLASctl")
if ("--benchmarks" %in% args) packages <- c(packages, "RGCCA", "PMA", "ggplot2")
if ("--checks" %in% args) packages <- c(packages, "testthat", "roxygen2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("Dependency check complete. The SGCA initializer is already bundled in egcar.\n")
