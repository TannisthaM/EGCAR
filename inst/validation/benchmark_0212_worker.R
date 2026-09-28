# Usage: Rscript benchmark_0212_worker.R library_directory output_prefix [reps]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
library(egcar)
reps <- if (length(args) >= 3L) as.integer(args[[3L]]) else 3L
set.seed(21205)
large <- length(args) >= 4L && args[[4L]] == "large"
n <- if (large) 30L else 50L
dims <- if (large) c(800L, 900L, 1000L) else c(320L, 360L, 400L)
latent <- matrix(rnorm(n * 2L), n, 2L)
X <- lapply(dims, function(p) {
  x <- matrix(rnorm(n * p), n, p)
  x[, 1:8] <- x[, 1:8] + latent %*% matrix(rnorm(16), 2L, 8L)
  x
})
ctl <- egcar_control(max_iter = if (large) 20L else 40L, abs_tol = 0, rel_tol = 0,
  check_every = 5L, adaptive_mu = TRUE, partial_eigen = TRUE)
p <- egcar_prepare(X)
rows <- list(); results <- list()
for (penalty in c("l11", "l21")) {
  for (i in seq_len(reps)) {
    # Each repetition includes construction of the loading factors.
    p$prep$loading_factor_cache <- new.env(parent = emptyenv())
    gc()
    elapsed <- system.time(fit <- egcar_fit(p, 2L, penalty, .02, ctl))[["elapsed"]]
    state_fields <- if (penalty == "l11") c("C", "Z", "H") else
      if (!is.null(fit$solver[["a"]])) c("C", "Hk", "Hl", "a") else c("C", "G", "V")
    rows[[length(rows) + 1L]] <- data.frame(version = as.character(packageVersion("egcar")),
      n = n, p = sum(dims), rank = 2L, penalty = penalty, rep = i,
      iterations = fit$iterations, fit_seconds = fit$fit_time,
      loading_seconds = fit$loading_time, total_seconds = elapsed,
      warm_state_bytes = as.numeric(object.size(fit$solver[state_fields])),
      loading_method = if (is.null(fit$loading$eigen_method)) "dense assembled operator with partial solve" else fit$loading$eigen_method)
    if (i == 1L) results[[penalty]] <- list(C = fit$C, U = fit$loading$U,
      values = fit$loading$eigenvalues, selected = fit$selected,
      primal = fit$solver$primal_residual, dual = fit$solver$dual_residual)
    rm(fit)
  }
}
write.csv(do.call(rbind, rows), paste0(args[[2L]], "_timings.csv"), row.names = FALSE)
saveRDS(results, paste0(args[[2L]], "_results.rds"))
writeLines(capture.output(sessionInfo()), paste0(args[[2L]], "_session.txt"))
