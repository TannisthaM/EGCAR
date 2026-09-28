#!/usr/bin/env Rscript
# Rscript tools/benchmark_0211.R /path/to/original/egcar /path/to/revised/egcar /path/to/results [reps]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Supply original source directory, revised source directory, output directory, and optional repetitions.")
sources <- normalizePath(args[1:2], mustWork = TRUE)
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[3L]], mustWork = TRUE)
reps <- if (length(args) > 3L) as.integer(args[[4L]]) else 3L
stopifnot(length(reps) == 1L, is.finite(reps), reps > 0L)
R <- file.path(R.home("bin"), "R")
Rscript <- file.path(R.home("bin"), "Rscript")
worker <- file.path(sources[[2L]], "inst", "validation", "benchmark_0211_worker.R")
labels <- c("original", "optimized")
for (i in seq_along(sources)) {
  lib <- file.path(out, paste0("lib_", labels[[i]])); dir.create(lib, showWarnings = FALSE)
  log <- file.path(out, paste0(labels[[i]], "_install.log"))
  status <- system2(R, c("CMD", "INSTALL", "--preclean", "--clean", "--no-multiarch",
    paste0("--library=", shQuote(lib)), shQuote(sources[[i]])), stdout = log, stderr = log)
  if (status != 0L) stop("Installation failed; inspect ", log)
}
# Alternate version order between repetitions to reduce order/thermal bias.
for (rep in seq_len(reps)) for (i in if (rep %% 2L) 1:2 else 2:1) {
  log <- file.path(out, paste0(labels[[i]], "_", rep, "_run.log"))
  prefix <- file.path(out, paste0(labels[[i]], "_", rep))
  status <- system2(Rscript, c(shQuote(worker),
    shQuote(file.path(out, paste0("lib_", labels[[i]]))), labels[[i]], shQuote(prefix), "1"),
    env = c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1",
            "VECLIB_MAXIMUM_THREADS=1"), stdout = log, stderr = log)
  if (status != 0L) stop("Benchmark failed; inspect ", log)
}
raw <- do.call(rbind, lapply(seq_len(reps), function(rep) do.call(rbind, lapply(labels, function(label) {
  z <- read.csv(file.path(out, paste0(label, "_", rep, "_timings.csv")))
  z$rep <- rep; z
}))))
write.csv(raw, file.path(out, "timings_all.csv"), row.names = FALSE)
med <- aggregate(cbind(elapsed_seconds, returned_object_bytes, serialized_bytes) ~ version + case + operation,
                 data = raw, FUN = median)
comparison <- merge(subset(med, version == "original", select = -version),
                    subset(med, version == "optimized", select = -version),
                    by = c("case", "operation"), suffixes = c("_original", "_optimized"))
comparison$speedup <- with(comparison, elapsed_seconds_original / elapsed_seconds_optimized)
comparison$returned_object_reduction_pct <- with(comparison,
  100 * (1 - returned_object_bytes_optimized / returned_object_bytes_original))
comparison$serialized_reduction_pct <- with(comparison,
  100 * (1 - serialized_bytes_optimized / serialized_bytes_original))
write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
a <- readRDS(file.path(out, "original_1_results.rds"))
b <- readRDS(file.path(out, "optimized_1_results.rds"))
close <- function(x, y, tol = 1e-6) isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
checks <- lapply(names(a), function(key) {
  x <- a[[key]]; y <- b[[key]]
  loading_distance <- if (is.null(x$L) && is.null(y$L)) 0 else if (is.null(x$L) || is.null(y$L)) Inf else {
    # Projection comparison avoids sign/rotation ambiguity of eigenvectors.
    qx <- qr.Q(qr(x$L)); qy <- qr.Q(qr(y$L))
    norm(tcrossprod(qx) - tcrossprod(qy), "F")
  }
  data.frame(case = key, C_equal = close(x$C, y$C), scores_equal = close(x$scores, y$scores),
    lambda_equal = close(x$lambda, y$lambda), loading_projector_distance = loading_distance,
    iterations_equal = identical(x$iterations, y$iterations),
    convergence_equal = identical(x$converged, y$converged), folds_equal = identical(x$fold_id, y$fold_id))
})
checks <- do.call(rbind, checks)
write.csv(checks, file.path(out, "equivalence.csv"), row.names = FALSE)
print(comparison, row.names = FALSE)
print(checks, row.names = FALSE)
if (!all(checks$C_equal & checks$scores_equal & checks$lambda_equal &
         checks$iterations_equal & checks$convergence_equal & checks$folds_equal &
         checks$loading_projector_distance < 1e-4)) stop("A numerical comparison failed; inspect equivalence.csv.")
cat("End-to-end comparison complete. Object/serialized sizes are not peak process RAM.\n")
