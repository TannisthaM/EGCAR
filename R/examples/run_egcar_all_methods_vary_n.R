#!/usr/bin/env Rscript

# =============================================================================
# EGCAR Midway3 vary_n -- built entirely on the installed egcar package's own
# exported functions: egcar::egcar_experiment_config() and
# egcar::run_egcar_experiments(). These are the SAME two functions used in
# run_EGCAR_package_signal08_r125.R. No estimator/CV/benchmark statistics are
# reimplemented in this file -- everything below this point is thin SLURM-
# array harness (config table + check_packages/install_packages/
# expected_tasks/task_info/worker/aggregate dispatch) needed to match
# run_egcar_all_methods_vary_n_array.sh, which is UNCHANGED and needs no
# edits to work with this file.
#
# Grid: n = {200,250,300,500,1000,10000} at p1=p2=p3=100, and
#       n = {250,300,450,1000,5000,10000} at p1=p2=p3=150 (unchanged),
# crossed with rank = 1:10 (was 1,2,5) and signal in {0.3,0.5,0.8}
# (was one fixed 0.8), 10 repetitions.
#   2 panels x 10 ranks x 3 signals x 6 n-values x 10 reps = 3,600 tasks.
#
# IMPORTANT ASSUMPTIONS THIS SCRIPT MAKES ABOUT egcar, WHICH I HAVE NOT BEEN
# ABLE TO VERIFY FROM HERE (no visibility into the installed package):
#   1. egcar::run_egcar_experiments() accepts a config scoped to a SINGLE n,
#      SINGLE rank, and SINGLE signal (n_grid/rank_grid of length 1) and
#      n_reps=1, and runs that one datapoint correctly rather than requiring
#      >=2 grid points.
#   2. Because n_reps=1 is used per call, this script gives each of the 10
#      repetitions of the SAME (n,rank,signal) its own master_seed (see
#      run_worker() below) so they are not identical runs. Whether this
#      actually produces independent draws depends on how egcar seeds
#      internally, which I cannot see.
#   3. Its output directory always contains a file literally named
#      simulation_results.csv with at least rank, n, method columns (as
#      required by replot_egcar_error_time() in the driver script you
#      uploaded).
#   4. It implements 10 methods (Oracle1-population, Oracle2-support,
#      EGCAR-L11-rate, EGCAR-L11-CV, EGCAR-L21-rate, EGCAR-L21-CV, SGCA,
#      RGCCA, SGCCA, MultiCCA) -- egcar_experiment_config()'s own argument
#      list has no "tied" penalty option, so that method (present in the
#      hand-rolled reference-engine version of this script) is NOT expected
#      here. If your egcar build has grown extra methods, they just appear
#      as extra rows; nothing below hardcodes exactly 10.
#
# RUN THE ONE-TASK SMOKE TEST BEFORE TRUSTING ANY OF THIS AT SCALE. If
# assumption 1 is wrong, it will fail loudly and immediately on that single
# task rather than burning 3,600 tasks' worth of compute first.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
`%||%` <- function(x, y) { if (is.null(x) || length(x) == 0L) y else x }
ARGS <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(ARGS)) ARGS[[1L]] else "help"
EXPERIMENT <- "vary_n"

R_LIB_USER <- path.expand(Sys.getenv("R_LIBS_USER", "~/Rlibs"))
dir.create(R_LIB_USER, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(R_LIB_USER, .libPaths())))

# =============================================================================
# 1. Package install/check -- egcar only, matching check_egcar_environment()
#    in run_EGCAR_package_signal08_r125.R (same minimum version, same
#    required CRAN packages), extended with an install_packages mode for
#    unattended SLURM use.
# =============================================================================

CRAN_PACKAGES <- c("future", "future.apply", "RGCCA", "PMA", "ggplot2", "remotes")
EGCAR_MIN_VERSION <- "0.2.5"

install_required_packages <- function() {
  missing <- CRAN_PACKAGES[!vapply(
    CRAN_PACKAGES, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) {
    message("Installing missing CRAN packages into ", R_LIB_USER, ": ",
            paste(missing, collapse = ", "))
    install.packages(
      missing, repos = "https://cloud.r-project.org", lib = R_LIB_USER,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = {
        requested <- suppressWarnings(as.integer(Sys.getenv("EGCAR_INSTALL_CORES", "2")))
        if (!is.finite(requested) || requested < 1L) requested <- 2L
        min(requested, 8L)
      }
    )
  }

  needs_install <- !requireNamespace("egcar", quietly = TRUE) ||
    identical(Sys.getenv("EGCAR_REINSTALL", "0"), "1") ||
    (requireNamespace("egcar", quietly = TRUE) &&
       utils::packageVersion("egcar") < EGCAR_MIN_VERSION)
  if (needs_install) {
    if ("egcar" %in% loadedNamespaces()) unloadNamespace("egcar")
    local_source <- Sys.getenv("EGCAR_LOCAL_SOURCE", "")
    github_repo <- Sys.getenv("EGCAR_GITHUB_REPO", "")

    if (nzchar(local_source) && (file.exists(local_source) || dir.exists(local_source))) {
      message("Installing egcar from local source: ", local_source)
      if (dir.exists(local_source)) {
        remotes::install_local(local_source, lib = R_LIB_USER, dependencies = NA, upgrade = "never")
      } else if (grepl("\\.zip$", local_source, ignore.case = TRUE)) {
        tmp <- tempfile("egcar_install_")
        dir.create(tmp)
        utils::unzip(local_source, exdir = tmp)
        pkg_dir <- file.path(tmp, "egcar")
        stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))
        install.packages(pkg_dir, repos = NULL, type = "source", lib = R_LIB_USER,
                          INSTALL_opts = c("--preclean", "--clean"))
      } else {
        install.packages(local_source, repos = NULL, type = "source", lib = R_LIB_USER)
      }
    } else if (nzchar(github_repo)) {
      message("Installing egcar from GitHub repository: ", github_repo)
      remotes::install_github(github_repo, lib = R_LIB_USER, dependencies = NA, upgrade = "never")
    } else {
      stop(
        "egcar is not installed (or is older than ", EGCAR_MIN_VERSION, "). ",
        "Set EGCAR_LOCAL_SOURCE to your zip/tarball/directory (the same ",
        "egcar_0.2.5_config_compat-style file you installed by hand), or ",
        "EGCAR_GITHUB_REPO to owner/repository, then rerun install_packages."
      )
    }
  }
  invisible(check_required_packages(stop_on_missing = TRUE))
}

check_required_packages <- function(stop_on_missing = TRUE) {
  packages <- c(CRAN_PACKAGES, "egcar")
  rows <- lapply(packages, function(pkg) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    data.frame(
      package = pkg, installed = ok,
      version = if (ok) as.character(utils::packageVersion(pkg)) else NA_character_,
      library = if (ok) find.package(pkg) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)

  missing <- tab$package[!tab$installed]
  version_ok <- TRUE
  if (requireNamespace("egcar", quietly = TRUE)) {
    version_ok <- utils::packageVersion("egcar") >= EGCAR_MIN_VERSION
    if (!version_ok) message("Installed egcar is older than ", EGCAR_MIN_VERSION, "; run install_packages.")
    needed_fns <- c("egcar_experiment_config", "run_egcar_experiments")
    missing_fns <- needed_fns[
      !vapply(needed_fns, exists, logical(1), envir = asNamespace("egcar"), inherits = FALSE)
    ]
    if (length(missing_fns)) {
      message("egcar is missing expected exported function(s): ", paste(missing_fns, collapse = ", "))
      message("Available egcar namespace symbols (first 40): ",
              paste(utils::head(sort(ls(asNamespace("egcar"))), 40L), collapse = ", "))
    }
  }
  if (stop_on_missing && (length(missing) || !version_ok)) {
    stop("Package preflight failed. Run this script with mode install_packages, then rerun check_packages.")
  }
  invisible(tab)
}

# =============================================================================
# 2. Config table -- same n/p_per_block/panel values as before; rank and
#    signal are the only dimensions that changed (1:10 and 3 levels).
# =============================================================================

N_REPS <- 10L
MASTER_SEED <- 20260907L

build_configs <- function() {
  A <- expand.grid(
    rank = 1:10,
    signal = c(0.3, 0.5, 0.8),
    n = c(200L, 250L, 300L, 500L, 1000L, 10000L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  A$p_per_block <- 100L
  A$panel <- "p1=p2=p3=100"

  B <- expand.grid(
    rank = 1:10,
    signal = c(0.3, 0.5, 0.8),
    n = c(250L, 300L, 450L, 1000L, 5000L, 10000L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  B$p_per_block <- 150L
  B$panel <- "p1=p2=p3=150"

  out <- rbind(A, B)
  out <- out[order(out$p_per_block, out$signal, out$rank, out$n),
             c("panel", "p_per_block", "n", "rank", "signal")]
  rownames(out) <- NULL
  out
}

CONFIGS <- build_configs()
CONFIGS$config_id <- seq_len(nrow(CONFIGS))
EXPECTED_TASKS <- nrow(CONFIGS) * N_REPS
TASKS <- CONFIGS[rep(seq_len(nrow(CONFIGS)), each = N_REPS), , drop = FALSE]
TASKS$rep_id <- rep(seq_len(N_REPS), times = nrow(CONFIGS))
TASKS$global_task <- seq_len(nrow(TASKS))
TASKS <- TASKS[, c("global_task", "config_id", "rep_id", "panel", "p_per_block", "n", "rank", "signal")]
rownames(TASKS) <- NULL

METHOD_ORDER <- c(
  "Oracle1-population", "Oracle2-support", "EGCAR-L11-rate", "EGCAR-L11-CV",
  "EGCAR-L21-rate", "EGCAR-L21-CV", "SGCA", "RGCCA", "SGCCA", "MultiCCA"
)
METHOD_COLORS <- setNames(
  grDevices::hcl.colors(length(METHOD_ORDER), palette = "Dark 3"), METHOD_ORDER
)

# =============================================================================
# 3. Small I/O helpers (unchanged in spirit from the reference-engine
#    version -- atomic writes, output directory layout).
# =============================================================================

atomic_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("atomic_", tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(x, tmp, row.names = row.names)
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
  invisible(path)
}

make_output_dirs <- function(outdir) {
  for (d in c("metrics", "raw", "completed", "plots_all_methods", "plots_without_oracle1"))
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# 4. worker(): the ONLY place egcar's functions are actually called. Builds
#    a config scoped to this task's single (n, rank, signal), requests
#    n_reps=1, and relabels egcar's own simulation_results.csv with the
#    bookkeeping columns (config_id, rep, panel, p_per_block, signal) the
#    aggregate step needs. No statistics are computed here -- that is all
#    inside egcar::run_egcar_experiments().
# =============================================================================

run_worker <- function(config_id, rep_id, outdir) {
  make_output_dirs(outdir)
  row <- CONFIGS[CONFIGS$config_id == config_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown config_id: ", config_id)
  if (rep_id < 1L || rep_id > N_REPS) stop("rep_id must be between 1 and ", N_REPS)

  tag <- sprintf("%d_%d", config_id, rep_id)
  done_path <- file.path(outdir, "completed", paste0("task_", tag, ".done"))
  if (file.exists(done_path) && !identical(Sys.getenv("EGCAR_OVERWRITE", "0"), "1")) {
    cat("Already complete; not recomputing: ", tag, "\n", sep = "")
    return(invisible(NULL))
  }

  panel_index <- match(row$panel[[1L]], sort(unique(CONFIGS$panel)))
  signal_index <- match(row$signal[[1L]], sort(unique(CONFIGS$signal)))
  # Distinct seed per (panel, signal, rank, rep) -- see assumption 2 in the
  # header: this is how the 10 repetitions of the same (n,rank,signal) get
  # different random draws, since each call below only asks egcar for one
  # repetition (n_reps=1), bypassing whatever internal rep-indexing egcar
  # would otherwise use for a multi-rep call.
  task_seed <- MASTER_SEED +
    (panel_index - 1L) * 100000000L +
    (signal_index - 1L) * 10000000L +
    rep_id * 100000L +
    row$rank[[1L]] * 1000L

  p_per_block <- row$p_per_block[[1L]]

  cfg <- egcar::egcar_experiment_config(
    p_list = rep(p_per_block, 3L),
    n_grid = row$n[[1L]],
    rank_grid = row$rank[[1L]],

    active_per_view = 5L,
    toeplitz_rho = c(0.5, 0.7, 0.9),
    signal = row$signal[[1L]],
    master_seed = task_seed,

    rho_e_cv_grid = 10^seq(-5, 4),
    lambda_g_cv_grid = 10^seq(-5, 4),
    rate_c_e = 1,
    rate_c_g = 1,
    n_folds = 5L,

    mu_z = 1,
    mu_g = 1,
    max_iter_cv = 1000L,
    max_iter_final = 2000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    adaptive_mu = TRUE,
    row_threshold = 1e-4,
    covariance_ridge = 1e-4,
    group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10,
    check_every_admm = 5L,

    run_external_benchmarks = TRUE,
    stop_if_benchmark_packages_missing = TRUE,

    sgca_k_grid = sort(unique(c(5L, 10L, 15L, 20L, 30L, 3L * p_per_block))),
    sgca_rho_grid = c(0, 1e-3, 1e-2, 0.1, 0.5, 1),
    sgca_lambda_grid = 10^seq(-5, 4),
    sgca_eta = 0.001,
    sgca_ridge_b = 1e-6,
    sgca_init_tol = 5e-3,
    sgca_max_iter_init = 1000L,
    sgca_tgd_tol = 1e-6,
    sgca_max_iter_tgd = 15000L,

    rgcca_tau_grid = c(1e-6, 1e-3, 0.1, 0.25, 0.5, 0.75, 1),
    rgcca_scheme = "factorial",
    rgcca_tol = 1e-8,
    rgcca_max_iter = 1000L,

    sgcca_sparsity_grid = seq(1 / sqrt(p_per_block), 1, length.out = 10L),
    multicca_l1_grid = seq(1, sqrt(p_per_block), length.out = 10L),
    multicca_niter = 25L,
    align_external_block_signs = TRUE,

    oracle1_max_iter = 2000L,
    fast_sgca_initializer = TRUE,
    multicca_backend = "gram",

    save_fits = FALSE,
    save_cv_fold_results = FALSE,
    save_loading_data = FALSE,
    save_compact_loadings = FALSE,
    retain_benchmark_fits = FALSE,
    loading_factor_cache_max = 4L,

    # This task's own plot would be a single point; aggregate_results()
    # below makes the real plots from every task's combined output.
    make_plots = FALSE,
    make_loading_plots = FALSE,
    loading_methods_per_page = 3L
  )

  requested_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6")))
  if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 6L
  workers <- max(1L, min(requested_cores - 1L, cfg$n_folds %||% 5L))

  raw_dir <- file.path(outdir, "raw", paste0("task_", tag))
  unlink(raw_dir, recursive = TRUE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf(
    "Starting %s: config=%d rep=%d n=%d p_per_block=%d rank=%d signal=%.1f; %d CV workers\n",
    EXPERIMENT, config_id, rep_id, row$n[[1L]], p_per_block, row$rank[[1L]], row$signal[[1L]], workers
  ))

  egcar::run_egcar_experiments(
    output_dir = raw_dir,
    workers = workers,
    n_reps = 1L,
    config = cfg,
    backend = "cpp",
    smoke_test = FALSE
  )

  results_path <- file.path(raw_dir, "simulation_results.csv")
  if (!file.exists(results_path)) {
    stop(
      "egcar::run_egcar_experiments() did not write the expected ",
      "simulation_results.csv at ", results_path, " -- check what file(s) ",
      "this installed version of egcar actually wrote: ",
      paste(list.files(raw_dir, recursive = TRUE), collapse = ", ")
    )
  }
  metrics <- utils::read.csv(results_path, stringsAsFactors = FALSE)
  metrics$config_id <- config_id
  metrics$rep <- rep_id
  metrics$panel <- row$panel[[1L]]
  metrics$p_per_block <- p_per_block
  metrics$p_total <- 3L * p_per_block
  metrics$signal <- row$signal[[1L]]
  metrics$experiment <- EXPERIMENT

  atomic_csv(metrics, file.path(outdir, "metrics", paste0("metrics_", tag, ".csv")))

  for (extra in c("cv_grid_results.csv", "selected_cv_parameters.csv")) {
    src <- file.path(raw_dir, extra)
    if (file.exists(src)) {
      d <- utils::read.csv(src, stringsAsFactors = FALSE)
      d$config_id <- config_id
      d$rep <- rep_id
      atomic_csv(d, file.path(outdir, "raw", paste0(tools::file_path_sans_ext(extra), "_", tag, ".csv")))
    }
  }

  writeLines(paste("completed", Sys.time()), done_path)
  cat("Completed: ", tag, "\n", sep = "")
  invisible(metrics)
}

# =============================================================================
# 5. aggregate_results(): combines every task's metrics_*.csv, writes the
#    summary table, and makes the (rank x signal) grid of plots -- all
#    methods and Oracle1-omitted, matching the plot naming/style used
#    throughout this project (fixed color per method, solid lines only).
# =============================================================================

aggregate_results <- function(outdir) {
  make_output_dirs(outdir)
  files <- list.files(file.path(outdir, "metrics"), pattern = "^metrics_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No metric files found in ", file.path(outdir, "metrics"))
  M <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  atomic_csv(M, file.path(outdir, "all_metrics.csv"))

  numeric_metrics <- intersect(
    c("C_relative_error", "subspace_euclidean", "subspace_sigma0",
      "total_time", "support_recall", "support_fdp"),
    names(M)
  )
  groups <- split(M, interaction(M$panel, M$p_per_block, M$n, M$rank, M$signal, M$method, drop = TRUE))
  S <- do.call(rbind, lapply(groups, function(g) {
    row <- g[1L, c("panel", "p_per_block", "p_total", "n", "rank", "signal", "method"), drop = FALSE]
    for (metric in numeric_metrics) {
      x <- g[[metric]][is.finite(g[[metric]])]
      row[[paste0(metric, "_mean")]] <- if (length(x)) mean(x) else NA_real_
      row[[paste0(metric, "_se")]] <- if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else
        if (length(x)) 0 else NA_real_
      row[[paste0(metric, "_n")]] <- length(x)
    }
    row
  }))
  atomic_csv(S, file.path(outdir, sprintf("summary_mean_over_%d_repetitions.csv", N_REPS)))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    make_plot_set <- function(results, omit_oracle1) {
      subdir <- if (omit_oracle1) "plots_without_oracle1" else "plots_all_methods"
      d <- if (omit_oracle1) results[results$method != "Oracle1-population", , drop = FALSE] else results
      output_dir <- file.path(outdir, subdir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      spec <- data.frame(
        metric = c("C_relative_error", "subspace_euclidean", "subspace_sigma0",
                   "total_time", "support_recall", "support_fdp"),
        label = c("Relative operator error", "Euclidean sine-theta distance",
                  "Sigma0 sine-theta distance", "Total time (seconds)",
                  "Row-support recall", "Row-support FDP"),
        log_y = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
        stringsAsFactors = FALSE
      )
      for (rank_value in sort(unique(d$rank)))
        for (signal_value in sort(unique(d$signal)))
          for (j in seq_len(nrow(spec))) {
            metric <- spec$metric[[j]]
            dd <- d[d$rank == rank_value & d$signal == signal_value & is.finite(d[[metric]]), , drop = FALSE]
            if (!nrow(dd)) next
            groups2 <- split(dd, interaction(dd$panel, dd$p_per_block, dd$n, dd$method, drop = TRUE))
            pd <- do.call(rbind, lapply(groups2, function(g) data.frame(
              panel = g$panel[[1L]], n = g$n[[1L]], method = as.character(g$method[[1L]]),
              mean = mean(g[[metric]]), stringsAsFactors = FALSE
            )))
            pd$method <- factor(pd$method, levels = METHOD_ORDER)
            pd$plot_value <- if (spec$log_y[[j]]) pmax(pd$mean, 1e-12) else pd$mean
            z <- ggplot2::ggplot(pd, ggplot2::aes(x = n, y = plot_value, color = method, group = method)) +
              ggplot2::geom_line(linewidth = 0.7, linetype = "solid", na.rm = TRUE) +
              ggplot2::geom_point(size = 2, shape = 16, na.rm = TRUE) +
              ggplot2::scale_x_log10() +
              ggplot2::scale_color_manual(values = METHOD_COLORS) +
              ggplot2::facet_wrap(~ panel, scales = "free_x", nrow = 1L) +
              ggplot2::labs(
                x = "Sample size n", y = spec$label[[j]],
                title = paste0(spec$label[[j]], ", r = ", rank_value, ", signal = ", signal_value,
                               if (omit_oracle1) " (Oracle1 omitted)" else " (all methods)"),
                color = "Method"
              ) +
              ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "bottom")
            if (spec$log_y[[j]]) z <- z + ggplot2::scale_y_log10()
            stem <- paste0(metric, "_rank_", rank_value, "_signal_", format(signal_value, nsmall = 1))
            ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), z, width = 12, height = 6.5, device = "pdf")
            ggplot2::ggsave(file.path(output_dir, paste0(stem, ".png")), z, width = 12, height = 6.5, dpi = 180)
          }
    }
    make_plot_set(M, FALSE)
    make_plot_set(M, TRUE)
  } else {
    message("ggplot2 not available; skipped plotting, all_metrics.csv/summary CSV were still written.")
  }

  writeLines(
    c(paste("Experiment:", EXPERIMENT),
      paste("Metric files found:", length(files), "of expected", EXPECTED_TASKS),
      paste("Generated:", Sys.time())),
    file.path(outdir, if (length(files) == EXPECTED_TASKS) "AGGREGATION_COMPLETE.txt" else "AGGREGATION_PARTIAL.txt")
  )
  if (length(files) < EXPECTED_TASKS) warning("Partial aggregation: ", length(files), "/", EXPECTED_TASKS, " tasks.")
  cat("Aggregation written to ", outdir, "\n", sep = "")
  invisible(M)
}

# =============================================================================
# 6. Dispatch -- same CLI shape run_egcar_all_methods_vary_n_array.sh already
#    expects: check_packages, install_packages, expected_tasks, task_info,
#    worker, aggregate.
# =============================================================================

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript SCRIPT.R check_packages\n",
    "  Rscript SCRIPT.R install_packages\n",
    "  Rscript SCRIPT.R expected_tasks\n",
    "  Rscript SCRIPT.R task_info <global_task>\n",
    "  Rscript SCRIPT.R worker <config_id> <rep_id> <out_dir>\n",
    "  Rscript SCRIPT.R aggregate 0 0 <out_dir>\n",
    sep = ""
  )
}

if (identical(MODE, "install_packages")) {
  install_required_packages()
} else if (identical(MODE, "check_packages")) {
  check_required_packages(TRUE)
} else if (identical(MODE, "expected_tasks")) {
  cat(EXPECTED_TASKS, "\n", sep = "")
} else if (identical(MODE, "task_info")) {
  if (length(ARGS) < 2L) stop("Missing global_task")
  id <- suppressWarnings(as.integer(ARGS[[2L]]))
  if (is.na(id) || id < 1L) stop("global_task must be a positive integer")
  z <- TASKS[TASKS$global_task == id, , drop = FALSE]
  if (nrow(z) != 1L) stop("global_task outside 1..", EXPECTED_TASKS)
  cat(z$config_id, z$rep_id, z$n, z$p_per_block, z$rank, z$signal, sep = "\t")
  cat("\n")
} else if (identical(MODE, "worker")) {
  if (length(ARGS) < 4L) stop("worker requires config_id, rep_id, out_dir")
  check_required_packages(TRUE)
  run_worker(as.integer(ARGS[[2L]]), as.integer(ARGS[[3L]]), ARGS[[4L]])
} else if (identical(MODE, "aggregate")) {
  if (length(ARGS) < 4L) stop("aggregate requires 0 0 out_dir")
  aggregate_results(ARGS[[4L]])
} else if (identical(MODE, "help")) {
  usage()
} else {
  usage(); stop("Unknown mode: ", MODE)
}
