#!/usr/bin/env Rscript
# Rscript tools/benchmark_0212.R old_R_library new_R_library output_dir [reps] [large]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Supply old/new R library directories and an output directory.")
libs <- vapply(args[1:2], normalizePath, character(1L), mustWork = TRUE)
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[3L]], mustWork = TRUE)
reps <- if (length(args) >= 4L) as.integer(args[[4L]]) else 3L
large <- length(args) >= 5L && args[[5L]] == "large"
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))))
worker <- file.path(root, "inst", "validation", "benchmark_0212_worker.R")
for (i in 1:2) {
  label <- c("old", "new")[[i]]
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(worker),
    shQuote(libs[[i]]), shQuote(file.path(out, label)), reps, if (large) "large"),
    stdout = file.path(out, paste0(label, ".log")), stderr = "")
  if (status != 0L) stop(label, " benchmark failed; inspect its log.")
}
a <- readRDS(file.path(out, "old_results.rds")); b <- readRDS(file.path(out, "new_results.rds"))
checks <- lapply(names(a), function(penalty) {
  x <- a[[penalty]]; y <- b[[penalty]]
  stopifnot(identical(x$selected, y$selected))
  # Rotations/signs within the requested subspace do not change the estimator.
  d <- data.frame(penalty = penalty, C_error = max(abs(unlist(x$C) - unlist(y$C))),
    eigenvalue_error = max(abs(x$values - y$values)),
    projector_error = max(abs(tcrossprod(x$U) - tcrossprod(y$U))),
    primal_error = abs(x$primal - y$primal), dual_error = abs(x$dual - y$dual))
  stopifnot(all(is.finite(as.matrix(d[-1L]))), all(as.matrix(d[-1L]) < 1e-7))
  d
})
write.csv(do.call(rbind, checks), file.path(out, "numerical_agreement.csv"), row.names = FALSE)
timings <- rbind(read.csv(file.path(out, "old_timings.csv")), read.csv(file.path(out, "new_timings.csv")))
summary <- aggregate(timings[c("fit_seconds", "loading_seconds", "total_seconds", "warm_state_bytes")],
  timings[c("version", "n", "p", "rank", "penalty")], median)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
