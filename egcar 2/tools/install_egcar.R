#!/usr/bin/env Rscript
# Rscript tools/install_egcar.R /path/to/egcar_0.2.0.tar.gz [--benchmarks]
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply a source archive or package-directory path.")
path <- normalizePath(args[[1L]], mustWork = TRUE)
packages <- c("Rcpp", "RcppArmadillo")
if ("--benchmarks" %in% args) packages <- c(packages,
  "future", "future.apply", "RGCCA", "PMA", "ggplot2", "RSpectra", "RhpcBLASctl")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
install.packages(path, repos = NULL, type = "source")
print(egcar::benchmark_dependencies(), row.names = FALSE)
