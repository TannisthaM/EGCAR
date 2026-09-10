# One small configuration, all ten methods, five shared folds / five CV workers.
# Install egcar >= 0.2.1 and its optional comparison/parallel dependencies first.
required <- c("egcar", "RcppEigen", "future", "future.apply", "RGCCA", "PMA")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
stopifnot(utils::packageVersion("egcar") >= "0.2.1")
library(egcar)

n <- 60L
p_total <- 12L                         # Three views, four variables per view.
r <- 1L
signal <- 0.8
out_dir <- "egcar_small_check"

# Small positive tuning grids; solver tolerances/iteration limits are unchanged.
cfg <- egcar_experiment_config(
  p_list = rep(p_total %/% 3L, 3L), n_grid = n, rank_grid = r,
  active_per_view = 2L, signal = signal, n_folds = 5L,
  master_seed = 20260910L,
  rho_e_cv_grid = c(0.001, 0.01), lambda_g_cv_grid = c(0.001, 0.01),
  sgca_k_grid = c(6L, 12L), sgca_rho_grid = c(0.01, 0.1),
  sgca_lambda_grid = c(0.1, 1), rgcca_tau_grid = c(0.5, 1),
  sgcca_sparsity_grid = c(0.75, 1), multicca_l1_grid = c(1.5, 2),
  run_external_benchmarks = TRUE,
  make_plots = FALSE, make_loading_plots = FALSE, save_fits = TRUE
)
stopifnot(p_total %% 3L == 0L, sum(cfg$p_list) == p_total)

ans <- run_egcar_experiments(
  output_dir = out_dir, workers = 5L, n_reps = 1L,
  config = cfg, backend = "cpp", smoke_test = FALSE
)
# smoke_test=FALSE prevents the runner from replacing this single configuration.

allocation <- read.csv(file.path(out_dir, "cv_worker_allocation.csv"))
a <- allocation[allocation$enabled, , drop = FALSE]
stopifnot(nrow(a) == 6L, all(a$allocated_cv_workers == 5L),
          all(a$folds == 5L), all(a$parallel_cv))

expected <- c("Oracle1-population", "Oracle2-support", "EGCAR-L11-rate",
  "EGCAR-L11-CV", "EGCAR-L21-rate", "EGCAR-L21-CV", "SGCA", "RGCCA", "SGCCA", "MultiCCA")
d <- ans$results
stopifnot(nrow(d) == length(expected), setequal(as.character(d$method), expected),
          all(d$n == n), all(d$rank == r))
report <- d[, c("method", "status", "converged", "subspace_sigma0", "total_time")]
report$finite_loading <- is.finite(d$subspace_euclidean) & is.finite(d$subspace_sigma0)
print(report, row.names = FALSE)
write.csv(report, file.path(out_dir, "small_check_report.csv"), row.names = FALSE)
if (!all(report$finite_loading)) {
  stop("Some methods did not return a usable loading. See small_check_report.csv, ",
       "egcar_cv_fold_results.csv and external_cv_fold_results.csv in ", out_dir, ".")
}
if (any(report$converged %in% FALSE)) {
  warning("All methods returned usable loadings, but some did not satisfy their stopping rules.")
}
message("All ten methods returned finite rank-one loadings. Convergence is reported separately.")
