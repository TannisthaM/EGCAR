#!/usr/bin/env Rscript
# Run the complete original local experiment through the installed egcar package.
# Terminal: Rscript run_EGCAR_from_package.R OUTPUT_DIR WORKERS N_REPS
# RStudio: source("run_EGCAR_from_package.R"); run_package_experiments(...)
# SGCA initialization is bundled; only RGCCA/PMA/ggplot2 are external comparisons/plots.
run_package_experiments <- function(
    output_dir = "egcar_package_outputs", workers = 1L, n_reps = 1L,
    smoke_test = identical(Sys.getenv("EGCAR_SMOKE_TEST", "0"), "1")) {
  if (!requireNamespace("egcar", quietly = TRUE)) stop("Install egcar first.")
  if (utils::packageVersion("egcar") < "0.2.0") stop("This runner requires egcar >= 0.2.0.")

  # These are the original full-study defaults, NOT a tiny illustrative example.
  cfg <- egcar::egcar_experiment_config(
    p_list = c(15L, 15L, 15L),
    n_grid = c(30L, 45L, 100L, 1000L, 5000L, 10000L),
    rank_grid = c(1L, 2L, 5L),
    active_per_view = 5L,
    toeplitz_rho = c(0.5, 0.7, 0.9),
    signal = 0.8,
    master_seed = 20260907L,
    n_folds = 5L,
    rho_e_cv_grid = 10^seq(-5, 4),
    lambda_g_cv_grid = 10^seq(-5, 4),
    run_external_benchmarks = TRUE,
    make_plots = TRUE,
    make_loading_plots = TRUE,
    save_fits = TRUE)

  # Additional controls are editable, for example:
  # cfg$max_iter_cv <- 1000L
  # cfg$max_iter_final <- 2000L
  # cfg$loading_plot_n <- c(100L, 10000L)
  # cfg$sgca_k_grid <- sort(unique(c(5L, 10L, 15L, 20L, 30L, sum(cfg$p_list))))

  egcar::run_egcar_experiments(output_dir = output_dir, workers = workers,
    n_reps = n_reps, config = cfg, backend = "cpp", smoke_test = smoke_test)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_package_experiments(
    output_dir = if (length(args) >= 1L) args[[1L]] else "egcar_package_outputs",
    workers = if (length(args) >= 2L) as.integer(args[[2L]]) else 1L,
    n_reps = if (length(args) >= 3L) as.integer(args[[3L]]) else 1L)
}
