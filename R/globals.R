# Namespace defaults and isolated per-call numerical environments.

EGCAR_BACKEND <- "cpp"
EGCAR_MASTER_BACKEND <- "cpp"
EGCAR_PARTIAL_EIGEN <- TRUE
EGCAR_PARTIAL_EIGEN_MIN <- 128L
LOADING_FACTOR_CACHE_MAX <- 4L
MU_Z <- 1
MU_G <- 1
ABS_TOL <- 1e-5
REL_TOL <- 1e-4
ADAPTIVE_MU <- TRUE
GROUP_ZERO_TOL <- 1e-8
ENTRY_ZERO_TOL <- 1e-10
ROW_THRESHOLD <- 1e-4
COVARIANCE_RIDGE <- 1e-4
ORACLE1_MAX_ITER <- 2000L
PARALLEL_CV <- FALSE
CV_WORKERS <- 1L
BLAS_THREADS <- 1L
SGCA_K_GRID <- c(5L, 10L, 15L, 20L, 30L, 45L)
SGCA_RHO_GRID <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1)
SGCA_LAMBDA_GRID <- 10^seq(-5, 4)
SGCA_ETA <- 0.001
SGCA_RIDGE_B <- 1e-6
SGCA_INIT_TOL <- 5e-3
SGCA_MAX_ITER_INIT <- 1000L
SGCA_TGD_TOL <- 1e-6
SGCA_MAX_ITER_TGD <- 15000L
FAST_SGCA_INITIALIZER <- TRUE
RGCCA_TAU_GRID <- c(1e-6, 1e-3, 0.1, 0.25, 0.5, 0.75, 1)
RGCCA_SCHEME <- "factorial"
RGCCA_TOL <- 1e-8
RGCCA_MAX_ITER <- 1000L
SGCCA_SPARSITY_GRID <- seq(1/sqrt(15), 1, length.out = 10L)
MULTICCA_L1_GRID <- seq(1, sqrt(15), length.out = 10L)
MULTICCA_NITER <- 25L
MULTICCA_BACKEND <- "gram"
ALIGN_EXTERNAL_BLOCK_SIGNS <- TRUE

.egcar_engine_functions <- c(
  "%||%",
  "symmetrize",
  "frob",
  "row_l2",
  "soft_threshold",
  "row_group_threshold",
  "matrix_power_psd",
  "block_diag",
  "make_block_indices",
  "edge_key",
  "make_edge_table",
  "make_incidence_layout",
  "empty_edge_list",
  "assemble_Mk",
  "assemble_all_M",
  "extract_edge_slice",
  "assemble_full_C",
  "split_full_C",
  "center_views",
  "center_views_at",
  "pairwise_loss",
  "operator_objective",
  "loading_metric_factors",
  "make_validation_covariance",
  "validation_score",
  "orthonormal_basis",
  "sine_theta_distance",
  "support_metrics",
  "make_folds",
  "make_fold_objects",
  "cache_problem_matrices_reference",
  "prepare_problem_reference",
  "fit_l11_admm_reference",
  "fit_l21_admm_reference",
  "fit_oracle_consensus_admm",
  "loading_from_operator",
  "egcar_project",
  "egcar_lift",
  "egcar_prepare_context",
  "cache_problem_matrices",
  "prepare_problem",
  "egcar_get_context",
  "egcar_initial_state",
  "egcar_view_copies",
  "egcar_group_norms",
  "egcar_loading_factors",
  "egcar_top_eigen",
  "egcar_loading_from_operator",
  "egcar_solve_R",
  "egcar_controls",
  "egcar_run_solver",
  "fit_l11_admm",
  "fit_l21_admm",
  "validation_loss",
  "validate_loading_matrix",
  "full_covariance_from_prep",
  "named_training_blocks",
  "training_gram_blocks",
  "synchronize_block_signs",
  "stack_block_weights",
  "bind_rows_fill",
  "summarize_loading_cv",
  "cross_validate_loading_grid",
  "sgca_hard_rows",
  "sgca_metric_normalize",
  "sgca_prepare_tgd_start",
  "sgca_prepare_initializer",
  "sgca_init_cached",
  "sgca_tgd_penalized",
  "get_sgca_initializer",
  "sgca_init_fixed",
  "updatePi",
  "updateH",
  "Soft",
  "sgca_common_cv",
  "rgcca_family_common_cv",
  "rgcca_common_cv",
  "sgcca_common_cv",
  "multicca_gram_fit",
  "multicca_common_cv",
  "toeplitz_correlation",
  "metric_normalize",
  "make_population",
  "rmvn_psd",
  "simulate_views",
  "prepare_population_problem",
  "fit_zero_penalty_oracle_admm",
  "fit_oracle_population",
  "fit_oracle_support"
)

set_blas_threads_one <- function() invisible(NULL)

parallel_map_candidates <- function(indices, FUN) lapply(indices, FUN)

.egcar_engine <- function(control = egcar_control(),
                          benchmarks = benchmark_control(), workers = 1L) {
  e <- new.env(parent = environment(.egcar_engine))
  for (nm in .egcar_engine_functions) {
    fn <- get(nm, envir = environment(.egcar_engine), inherits = FALSE)
    environment(fn) <- e
    assign(nm, fn, envir = e)
  }
  e$EGCAR_BACKEND <- tolower(control$backend)
  e$EGCAR_MASTER_BACKEND <- e$EGCAR_BACKEND
  e$EGCAR_PARTIAL_EIGEN <- control$partial_eigen
  e$EGCAR_PARTIAL_EIGEN_MIN <- control$partial_eigen_min
  e$LOADING_FACTOR_CACHE_MAX <- control$loading_cache_max
  e$MU_Z <- e$MU_G <- control$mu
  e$ABS_TOL <- control$abs_tol
  e$REL_TOL <- control$rel_tol
  e$ADAPTIVE_MU <- control$adaptive_mu
  e$GROUP_ZERO_TOL <- control$group_zero_tol
  e$ENTRY_ZERO_TOL <- control$entry_zero_tol
  e$ROW_THRESHOLD <- control$row_threshold
  e$COVARIANCE_RIDGE <- control$covariance_ridge
  e$ORACLE1_MAX_ITER <- control$max_iter
  e$CV_WORKERS <- as.integer(workers)
  e$PARALLEL_CV <- workers > 1L
  e$BLAS_THREADS <- control$blas_threads
  bm <- c(sgca_eta = "SGCA_ETA", sgca_ridge = "SGCA_RIDGE_B",
          sgca_init_tol = "SGCA_INIT_TOL", sgca_init_max_iter = "SGCA_MAX_ITER_INIT",
          sgca_tgd_tol = "SGCA_TGD_TOL", sgca_tgd_max_iter = "SGCA_MAX_ITER_TGD",
          fast_sgca_initializer = "FAST_SGCA_INITIALIZER",
          rgcca_scheme = "RGCCA_SCHEME", rgcca_tol = "RGCCA_TOL",
          rgcca_max_iter = "RGCCA_MAX_ITER", multicca_niter = "MULTICCA_NITER",
          multicca_backend = "MULTICCA_BACKEND",
          align_signs = "ALIGN_EXTERNAL_BLOCK_SIGNS")
  for (nm in names(bm)) assign(bm[[nm]], benchmarks[[nm]], envir = e)
  # Thread limits are scoped by the public entry points and dispatcher.
  e$set_blas_threads_one <- function() invisible(NULL)
  dispatcher <- function(indices, FUN) {
    task <- function(i) .egcar_with_threads(BLAS_THREADS, FUN(i))
    if (PARALLEL_CV) {
      future.apply::future_lapply(indices, task, future.seed = TRUE,
        future.scheduling = 1, future.packages = "egcar")
    } else lapply(indices, task)
  }
  environment(dispatcher) <- e
  e$parallel_map_candidates <- dispatcher
  e
}
