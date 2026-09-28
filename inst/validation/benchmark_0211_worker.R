# Invoked in separate R processes by tools/benchmark_0211.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
.libPaths(c(args[[1L]], .libPaths()))
library(egcar)
label <- args[[2L]]; output <- args[[3L]]; reps <- as.integer(args[[4L]])
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
}
rows <- list(); results <- list()
record <- function(case, operation, rep, expr) {
  gc(reset = TRUE)
  t <- system.time(value <- force(expr))[["elapsed"]]
  # Object size may count shared R matrices repeatedly. Serialized size is
  # also useful for worker-transfer/storage costs; neither is process peak RSS.
  rows[[length(rows) + 1L]] <<- data.frame(version = label, case = case,
    operation = operation, rep = rep, elapsed_seconds = t,
    returned_object_bytes = as.numeric(object.size(value)),
    serialized_bytes = length(serialize(value, NULL)))
  value
}
# Data generation avoids timing the simulation/population oracle machinery.
cases <- list(tall = c(n = 180L, p = 45L), wide = c(n = 50L, p = 160L))
for (case in names(cases)) {
  n <- cases[[case]][["n"]]; p <- cases[[case]][["p"]]
  set.seed(9211)
  H <- matrix(rnorm(n * 2L), n, 2L)
  X <- lapply(c(p, p + 3L, p - 2L), function(pk) {
    Z <- matrix(rnorm(n * pk), n, pk)
    Z[, 1:6] <- Z[, 1:6] + H %*% matrix(rnorm(12), 2L, 6L)
    Z
  })
  fold_id <- rep(1:3, length.out = n)
  ctl <- egcar_control(backend = "cpp", max_iter = 120L, max_iter_cv = 80L,
    abs_tol = 0, rel_tol = 0, partial_eigen = FALSE)
  # Equal fixed iteration counts isolate computation from stopping variation.
  invisible(egcar_fit(X, 2L, "l11", .02, ctl))
  for (rep in seq_len(reps)) {
    prepared <- record(case, "prepare", rep, egcar_prepare(X))
    folds <- record(case, "cv_prepare", rep, egcar_cv_data(X, fold_id))
    for (penalty in c("l11", "l21")) {
      fit <- record(case, paste0(penalty, "_fit_prepared"), rep,
        egcar_fit(prepared, 2L, penalty, .02, ctl))
      cv <- record(case, paste0(penalty, "_cv_prepared"), rep,
        egcar_cv(folds, 2L, penalty, c(.04, .02, .005), workers = 1L, control = ctl))
      if (rep == 1L) results[[paste(case, penalty)]] <- list(C = fit$C, L = fit$L,
        status = fit$status, iterations = fit$iterations, converged = fit$converged,
        scores = cv$cv_fold_table$score, lambda = cv$lambda, fold_id = cv$fold_id)
    }
  }
}
write.csv(do.call(rbind, rows), paste0(output, "_timings.csv"), row.names = FALSE)
saveRDS(results, paste0(output, "_results.rds"))
writeLines(c(capture.output(sessionInfo()), paste("BLAS:", extSoftVersion()[["BLAS"]])),
           paste0(output, "_session.txt"))
