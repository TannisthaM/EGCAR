# RGCCA comparison and shared RGCCA-family fitting machinery.

#' @rdname comparison_cv
#' @export
rgcca_cv <- function(x, rank = 1L, tau_grid = c(1e-6, 1e-3, 0.1, 0.25, 0.5, 0.75, 1),
                     fold_id = NULL, nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control(), benchmarks = benchmark_control()) {
  .egcar_external_cv(x, rank, "RGCCA", fold_id, nfolds, seed, workers, control,
                     benchmarks, list(tau_grid = tau_grid), match.call())
}

rgcca_family_common_cv <- function(
    full_views, full_prep, fold_objects, rank, parameter_grid,
    method = c("rgcca", "sgcca"), seed = 1L,
    parallel_folds = PARALLEL_CV) {

  method <- match.arg(method)
  K <- length(full_prep$p_list)
  if (!requireNamespace("RGCCA", quietly = TRUE)) stop("This benchmark requires RGCCA.")
  if (!length(parameter_grid) || any(!is.finite(parameter_grid))) stop("Invalid RGCCA-family grid.")
  if (method == "rgcca") {
    if (any(parameter_grid < 0 | parameter_grid > 1)) stop("RGCCA tau must be in [0, 1].")
    grid <- data.frame(rgcca_tau = sort(unique(parameter_grid), decreasing = TRUE))
  } else {
    lower <- max(1 / sqrt(full_prep$p_list))
    if (any(parameter_grid < lower - 1e-12 | parameter_grid > 1)) {
      stop("Each tied SGCCA sparsity coefficient must be between max_j(1/sqrt(p_j)) and 1.")
    }
    grid <- data.frame(sgcca_sparsity = sort(unique(pmax(lower, parameter_grid))))
  }
  # Ties: stronger shrinkage (larger tau), or stronger sparsity (smaller bound).
  prepare_context <- function(train_views, prep, final) {
    list(blocks = named_training_blocks(train_views), prep = prep,
         gram_blocks = training_gram_blocks(prep))
  }
  fit_candidate <- function(context, parameters, final) {
    args <- list(
      blocks = context$blocks, connection = 1 - diag(K),
      method = method, ncomp = rep(as.integer(rank), K),
      scheme = RGCCA_SCHEME,
      scale = FALSE, scale_block = FALSE, bias = TRUE,
      init = "svd", verbose = FALSE, quiet = TRUE,
      tol = RGCCA_TOL, n_iter_max = RGCCA_MAX_ITER
    )
    if (method == "rgcca") {
      args$tau <- rep(parameters$rgcca_tau[[1L]], K)
    } else {
      # SGCCA tunes sparsity, NOT tau. This is the actual RGCCA API argument.
      args$sparsity <- rep(parameters$sgcca_sparsity[[1L]], K)
    }
    fit <- do.call(RGCCA::rgcca, args)
    # astar (NOT a) maps the ORIGINAL centered blocks to all deflated scores.
    L <- stack_block_weights(fit$astar, context, rank)
    # Package crit histories provide a limit diagnostic, not a guaranteed
    # optimizer convergence certificate. Do not mark every returned fit TRUE.
    hist <- fit$crit
    lens <- if (is.list(hist)) lengths(hist) else length(hist)
    hit_limit <- length(lens) > 0L && any(lens >= RGCCA_MAX_ITER)
    list(L = L, fit = fit,
         converged = if (hit_limit) FALSE else NA,
         iterations = if (length(lens)) as.integer(sum(lens)) else NA_integer_)
  }
  cross_validate_loading_grid(
    full_views, full_prep, fold_objects, rank, grid,
    prepare_context, fit_candidate, label = toupper(method), seed = seed,
    parallel_folds = parallel_folds
  )
}

rgcca_common_cv <- function(
    full_views, full_prep, fold_objects, rank,
    tau_grid = RGCCA_TAU_GRID, seed = 1L, parallel_folds = PARALLEL_CV) {
  rgcca_family_common_cv(
    full_views, full_prep, fold_objects, rank, tau_grid,
    method = "rgcca", seed = seed, parallel_folds = parallel_folds
  )
}
