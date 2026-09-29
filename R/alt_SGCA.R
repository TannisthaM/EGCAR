# Paper-based SGCA: cached initialization, fixed-step TGD and selectable CV score.

#' Cross-validation for SGCA, RGCCA, SGCCA and MultiCCA
#'
#' @description Paper-based SGCA and optional comparison wrappers using shared folds and fold-level parallelism.
#' @param x Raw views or the same \code{egcar_cv_data} object used for EGCAR.
#' @param rank Target loading rank.
#' @param k_grid SGCA total row-support budget across all views. NULL uses valid members of 5, 10, ..., 100; p is used if none are valid.
#' @param rho_grid NULL uses the paper rate per training sample; numeric values are direct initializer coefficients.
#' @param lambda_grid SGCA strictly positive penalized-gradient coefficients.
#' @param tau_grid RGCCA tied tau grid in [0,1].
#' @param sparsity_grid SGCCA tied sparsity bounds between the maximum of 1/sqrt(p_k) and one. NULL uses a dimension-valid grid.
#' @param penalty_grid MultiCCA tied L1 bounds between one and the minimum of sqrt(p_k). NULL uses a dimension-valid grid.
#' @param fold_id,nfolds,seed,workers Common-fold and parallel controls as in \code{egcar_cv}.
#' @param control EGCAR infrastructure controls, including scoped BLAS threads. EGCAR ADMM controls do not change comparison solvers.
#' @param benchmarks A \code{benchmark_control} object or named argument list.
#' @details SGCA defaults follow Gao and Ma (2023), Algorithm 1 and Sections 5.1--5.2:
#' fixed 15000 thresholded-gradient updates, eta=0.001, lambda=0.01, a generalized
#' Fantope initializer with trace exactly rank, and rho=0.5*sqrt(log(p)/n_train).
#' Covariances use divisor n. Initial loading extraction uses leading algebraic
#' eigenvectors, followed by hard row thresholding, metric normalization and
#' rescaling. Each gradient update uses 2*eta and has no backtracking or intermediate
#' normalization. The final loading is metric-normalized.
#' 
#' The initializer uses the authors' R-code tolerance 0.005 and an enforced
#' 1000-outer-iteration cap; these are implementation controls, not numerical
#' accuracy guarantees specified by the paper. The original R counter bug is not
#' reproduced. Initialization is cached within each fold and rho value.
#' 
#' By default only sparsity is tuned, on valid members of 5,10,...,100, using the
#' paper's held-out covariance trace score. On dimensions below the first valid
#' grid point, the public NULL-grid default uses all p rows. NULL rho_grid applies
#' the rate separately on each training fold and full refit; explicit numeric
#' rho_grid values retain their old DIRECT-coefficient meaning. A supplied
#' lambda_grid also overrides the paper default. Five folds are the default;
#' precomputed folds are respected.
#' 
#' benchmark_control(sgca_cv_score="common_loss") selects EGCAR's common score
#' instead; this is an explicit comparison variant. sgca_stopping="absolute_change"
#' allows early stopping at an absolute Frobenius iterate change of 1e-6 by default.
#' That variant still uses the paper's normalization and 2*eta update, so it is not
#' a byte-for-byte reproduction of the authors' R implementation.
#' 
#' For fixed iterations, status="ok" and completed=TRUE mean a finite full-rank
#' benchmark output completed the requested steps. converged separately requires
#' both numerical change tests. A finite capped initializer is allowed by the
#' benchmark default and its cap is reported. Use sgca_cv_require_convergence=TRUE
#' to reject such candidates. EGCAR's cv_require_convergence does not control SGCA.
#' Fold diagnostics and final diagnostics record actual rho, stage iterations,
#' change magnitudes, convergence and stop reasons, even when dense fits are dropped.
#' 
#' The separate six-hour budget covers the entire SGCA CV grid and refit, including
#' initialization, not each fit. A timeout returns status="time_limit" and no loading.
#' Shared data preparation precedes the budget. Finite budgets use a supervised R
#' process; sgca_time_limit=Inf disables supervision. Numerical failures are not
#' repaired by changing rank, adding an unrequested ridge or switching algorithms.#'
#' RGCCA and SGCCA use the public \code{RGCCA::rgcca} solver with variable/block
#' scaling disabled. They use \code{astar} loadings mapping original centered data
#' to deflated components. MultiCCA uses the supplied cached-Gram implementation
#' and the installed PMA threshold helpers, or \code{PMA::MultiCCA} when selected.
#' External block signs are resolved using training observations only.
#'
#' RGCCA, SGCCA and MultiCCA use the common EGCAR validation loss. SGCA uses
#' its selected score. A candidate
#' must have finite loss on every fold. A returned finite approximate fit does
#' not imply convergence; diagnostic flags are retained. These methods return
#' loading estimates, not an EGCAR cross-view operator.
#'
#' Optional packages are never automatically installed. Run the included
#' \code{00_install_dependencies.R --benchmarks} example explicitly when needed.
#' For external methods loading extraction is included in \code{fit_time};
#' \code{loading_time} is therefore zero, not an additional omitted cost.
#' @return An \code{egcar_cv} object with \code{L}, \code{best}, the external \code{fit_full}, candidate and fold diagnostics, status, means, controls, labels and timings.
#' @rdname comparison_cv
#' @examples
#' \dontrun{
#' sim <- egcar_simulate(n = 60)
#' shared <- egcar_cv_data(sim$views, nfolds = 3)
#' rgcca_cv(shared, tau_grid = c(0.1, 1))
#' sgcca_cv(shared, sparsity_grid = c(0.5, 1))
#' multicca_cv(shared, penalty_grid = c(2, 3))
#' sgca_cv(shared)
#' }
#' @export
sgca_cv <- function(x, rank = 1L, k_grid = NULL,
                    rho_grid = NULL,
                    lambda_grid = 0.01, fold_id = NULL, nfolds = 5L,
                    seed = 1L, workers = 1L, control = egcar_control(),
                    benchmarks = benchmark_control()) {
  .egcar_external_cv(x, rank, "SGCA", fold_id, nfolds, seed, workers, control,
    benchmarks, list(k_grid = k_grid, rho_grid = rho_grid, lambda_grid = lambda_grid), match.call())
}

sgca_hard_rows <- function(U, k) {
  U <- as.matrix(U)
  if (k < ncol(U) || k > nrow(U)) stop("SGCA requires rank <= k <= p.")
  if (k < nrow(U)) {
    # Definition 6: break ties by the smaller row index, including rank > 1.
    keep <- order(-rowSums(U * U), seq_len(nrow(U)))[seq_len(k)]
    U[setdiff(seq_len(nrow(U)), keep), ] <- 0
  }
  U
}

sgca_metric_normalize <- function(U, B, tol = 1e-10) {
  G <- symmetrize(crossprod(U, B %*% U))
  ev <- eigen(G, symmetric = TRUE)
  d <- ev$values
  if (any(!is.finite(d)) || d[[1L]] <= 0 ||
      min(d) <= tol * d[[1L]]) stop("SGCA has a rank-deficient training Gram matrix.")
  tcrossprod(sweep(U %*% ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
}

sgca_prepare_tgd_start <- function(A, B, init, k) {
  U0 <- sgca_metric_normalize(sgca_hard_rows(init, k), B)
  ev <- eigen(symmetrize(crossprod(U0, A %*% U0)), symmetric = TRUE)
  list(vectors = ev$vectors, values = ev$values,
       U0_vectors = U0 %*% ev$vectors)
}

sgca_prepare_initializer <- function(A, B, reference, nu = 1) {
  A <- symmetrize(A)
  B <- symmetrize(B)
  h_update <- tryCatch(get("updateH", envir = environment(reference), inherits = TRUE),
                       error = function(e) NULL)
  if (!is.function(h_update)) return(NULL)  # keep compatibility with other versions
  ev <- eigen(B, symmetric = TRUE)
  d <- pmax(ev$values, 0)
  sqB <- tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors)
  tau <- 4 * nu * max(d)^2
  if (!is.finite(tau) || tau <= 0) tau <- 1
  list(p = nrow(B), B = B, sqB = sqB, tau = tau, nu = nu,
       A_scaled = (1 / tau) * A, B_scale = nu / tau, update_H = h_update)
}

sgca_init_cached <- function(prepared, rho, K,
                             epsilon = SGCA_INIT_TOL,
                             maxiter = SGCA_MAX_ITER_INIT, trace = FALSE) {
  if (is.null(prepared)) stop("Missing SGCA initializer preparation.")
  z <- prepared
  p <- z$p
  H <- Pi <- oldPi <- diag(1, p)
  Gamma <- matrix(0, p, p)
  criteria <- Inf
  iter <- 0L
  threshold <- rho / z$tau
  while (criteria > epsilon && iter < maxiter) {
    fixed_H <- z$B_scale * (z$sqB %*% (H - Gamma / z$nu) %*% z$sqB)
    for (j in seq_len(20L)) {
      Pi <- soft_threshold(Pi + z$A_scaled -
        z$B_scale * (z$B %*% Pi %*% z$B) + fixed_H, threshold)
    }
    H <- z$update_H(z$sqB, Gamma, z$nu, Pi, K)
    Gamma <- Gamma + (z$sqB %*% Pi %*% z$sqB - H) * z$nu
    criteria <- sqrt(sum((Pi - oldPi)^2))
    oldPi <- Pi
    iter <- iter + 1L
    if (trace) cat("iter:", iter, "crit:", criteria, "\n")
  }
  list(Pi = Pi, H = H, Gamma = Gamma, iteration = iter, convergence = criteria)
}

# Algorithm 1 of Gao and Ma (2023), in the paper's unscaled V coordinates.
# Fixed T steps are the default; the tolerance remains a separate diagnostic.
sgca_tgd_penalized <- function(
    A, B, init, rank, k, lambda,
    eta = SGCA_ETA, max_iter = SGCA_MAX_ITER_TGD,
    tol = SGCA_TGD_TOL, stopping = SGCA_STOPPING,
    prepared_start = NULL, matrices_prepared = FALSE) {
  .egcar_scalar(lambda, "lambda", strict = TRUE)
  .egcar_scalar(eta, "eta", strict = TRUE)
  .egcar_scalar(tol, "tol", strict = TRUE)
  .egcar_optional_limit(max_iter, "max_iter", integer = TRUE)
  stopping <- match.arg(stopping, c("fixed_iterations", "absolute_change"))
  if (stopping == "fixed_iterations" && !is.finite(max_iter))
    stop("Fixed-iteration SGCA requires a finite max_iter.")
  if (!isTRUE(matrices_prepared)) {
    A <- symmetrize(A)
    B <- symmetrize(B)
  }
  if (is.null(prepared_start)) prepared_start <- sgca_prepare_tgd_start(A, B, init, k)
  st <- prepared_start
  # Lines 1--2: normalize in B, then multiply by (I + U' A U / lambda)^1/2.
  scales <- 1 + st$values / lambda
  if (any(!is.finite(scales)) || any(scales <= 0))
    stop("SGCA rescaling is not positive definite.")
  V <- tcrossprod(sweep(st$U0_vectors, 2L, sqrt(scales), "*"), st$vectors)
  I <- diag(rank)
  change <- Inf
  iter <- 0
  while (iter < max_iter) {
    AV <- A %*% V
    BV <- B %*% V
    # The factor 2 is part of the gradient of equation (13).
    gradient_half <- -AV + lambda * BV %*% (crossprod(V, BV) - I)
    next_V <- sgca_hard_rows(V - (2 * eta) * gradient_half, k)
    if (any(!is.finite(next_V)))
      stop("Non-finite SGCA iterate under the fixed paper step; no backtracking was applied.")
    change <- frob(next_V - V)
    V <- next_V
    iter <- iter + 1
    if (stopping == "absolute_change" && change <= tol) break
  }
  # Final normalization is included in the surrounding fit timer.
  L <- validate_loading_matrix(sgca_metric_normalize(V, B), nrow(A), rank)
  converged <- is.finite(change) && change <= tol
  completed <- if (stopping == "fixed_iterations") iter == max_iter else converged
  value <- -sum(V * (A %*% V)) + lambda / 2 * sum((crossprod(V, B %*% V) - I)^2)
  list(L = L, converged = converged, completed = completed, iterations = iter,
       objective = value, final_step = 2 * eta, absolute_change = change,
       stopping = stopping,
       stop_reason = if (stopping == "fixed_iterations") "fixed_iterations" else
         if (converged) "absolute_change" else "iteration_limit",
       lambda = lambda, k = k)
}

# Section 5.2: validation covariance is centered within the held-out sample;
# it is not re-normalized by the held-out block metric. Only r score columns
# are centered, avoiding a second n_test by p copy or a p by p covariance.
sgca_paper_validation_loss <- function(L, validation) {
  if (is.null(validation$n) || validation$n < 2L)
    stop("Paper SGCA scoring needs at least two validation observations; recreate older CV-data objects with egcar 0.2.16.")
  if (is.null(validation$views)) {
    if (is.null(validation$mean)) stop("Recreate the CV-data object with egcar 0.2.16 for paper SGCA scoring.")
    return(-(sum(L * (validation$Sigma %*% L)) - sum(crossprod(validation$mean, L)^2)))
  }
  scores <- Reduce(`+`, lapply(seq_along(validation$views), function(k)
    validation$views[[k]] %*% L[validation$indices[[k]], , drop = FALSE]))
  scores <- sweep(scores, 2L, colMeans(scores), "-")
  -sum(scores * scores) / validation$n
}

get_sgca_initializer <- function() {
  # Bundled initializer and its dependency closure: no optional namespace lookup.
  sgca_init_fixed
}

sgca_common_cv_run <- function(
    full_views, full_prep, fold_objects, rank,
    k_grid = SGCA_K_GRID, rho_grid = SGCA_RHO_GRID,
    lambda_grid = SGCA_LAMBDA_GRID, seed = 1L,
    parallel_folds = PARALLEL_CV) {

  p <- full_prep$p
  if (any(!is.finite(k_grid)) || any(k_grid != floor(k_grid))) {
    stop("SGCA k_grid must contain finite integers.")
  }
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[k_grid >= rank & k_grid <= p]
  if (!length(k_grid)) stop("No SGCA k satisfies rank <= k <= p.")
  if (!is.null(rho_grid) && (!length(rho_grid) || any(!is.finite(rho_grid)) || any(rho_grid < 0))) {
    stop("SGCA rho_grid must contain nonnegative, finite DIRECT coefficients.")
  }
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    stop("SGCA lambda_grid must contain positive, finite coefficients.")
  }
  grid <- expand.grid(
    sgca_k = k_grid, sgca_rho = if (is.null(rho_grid)) NA_real_ else sort(unique(rho_grid), decreasing = TRUE),
    sgca_lambda = sort(unique(lambda_grid)), KEEP.OUT.ATTRS = FALSE
  )
  grid$sgca_rho_rule <- if (is.null(rho_grid)) "0.5*sqrt(log(p)/n_train)" else "direct"
  grid$sgca_cv_score <- SGCA_CV_SCORE
  # Exact CV ties: smaller k, then larger rho, then smaller lambda.
  grid <- grid[order(grid$sgca_k, -grid$sgca_rho, grid$sgca_lambda), , drop = FALSE]
  rownames(grid) <- NULL

  prepare_context <- function(train_views, prep, final) {
    A <- full_covariance_from_prep(prep)
    B <- prep$Sigma0 %||% block_diag(prep$S_kk)
    diag(B) <- diag(B) + SGCA_RIDGE_B
    initializer <- get_sgca_initializer()
    list(
      prep = prep, A = A, B = B, initializer = initializer,
      initializer_prepared = if (FAST_SGCA_INITIALIZER)
        sgca_prepare_initializer(A, B, initializer) else NULL,
      init_cache = new.env(parent = emptyenv()),
      tgd_start_cache = new.env(parent = emptyenv())
    )
  }

  fit_candidate <- function(context, parameters, final) {
    rho <- if (is.null(rho_grid)) 0.5 * sqrt(log(p) / context$prep$n) else parameters$sgca_rho[[1L]]
    k <- parameters$sgca_k[[1L]]
    lambda <- parameters$sgca_lambda[[1L]]
    key <- sprintf("rho_%.17g", rho)
    if (!exists(key, envir = context$init_cache, inherits = FALSE)) {
      initialized <- tryCatch({
        z <- if (!is.null(context$initializer_prepared)) {
          sgca_init_cached(context$initializer_prepared, rho = rho, K = rank,
            epsilon = SGCA_INIT_TOL, maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
        } else {
          context$initializer(A = context$A, B = context$B, rho = rho, K = rank,
            nu = 1, epsilon = SGCA_INIT_TOL,
            maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
        }
        Pi <- as.matrix(z$Pi)
        if (any(!is.finite(Pi))) stop("Non-finite SGCA initializer.")
        # Equation (15) uses largest algebraic eigenvalues, not singular values.
        ss <- eigen(symmetrize(Pi), symmetric = TRUE)
        if (ss$values[[1L]] <= 0 || ss$values[[rank]] <= 1e-10 * ss$values[[1L]])
          stop("The SGCA initializer has fewer than rank positive leading eigenvalues.")
        U <- sweep(ss$vectors[, seq_len(rank), drop = FALSE], 2L,
                   sqrt(ss$values[seq_len(rank)]), "*")
        list(U = U, convergence = z$convergence, iteration = z$iteration,
             raw = if (isTRUE(final) && isTRUE(get0("RETAIN_BENCHMARK_FITS", inherits = TRUE, ifnotfound = TRUE))) z else NULL,
             error = NULL)
      }, error = function(e) list(U = NULL, convergence = NA_real_, iteration = NA_integer_,
                                  raw = NULL, error = conditionMessage(e)))
      assign(key, initialized, envir = context$init_cache)
    }
    ini <- get(key, envir = context$init_cache, inherits = FALSE)
    if (is.null(ini$U)) stop(ini$error)
    start_key <- paste0(key, "_k", k)
    if (!exists(start_key, envir = context$tgd_start_cache, inherits = FALSE)) {
      st <- tryCatch(list(value = sgca_prepare_tgd_start(context$A, context$B, ini$U, k)),
                     error = function(e) list(error = conditionMessage(e)))
      assign(start_key, st, envir = context$tgd_start_cache)
    }
    st <- get(start_key, envir = context$tgd_start_cache, inherits = FALSE)
    if (is.null(st$value)) stop(st$error)
    tgd <- sgca_tgd_penalized(
      A = context$A, B = context$B, init = ini$U,
      rank = rank, k = k, lambda = lambda,
      prepared_start = st$value, matrices_prepared = TRUE
    )
    init_conv <- if (is.numeric(ini$convergence) && length(ini$convergence) == 1L) {
      is.finite(ini$convergence) && ini$convergence <= SGCA_INIT_TOL
    } else NA
    conv <- if (identical(init_conv, FALSE)) FALSE else tgd$converged
    init_iters <- as.integer(ini$iteration %||% NA_integer_)
    keep_fit <- isTRUE(final) && isTRUE(get0("RETAIN_BENCHMARK_FITS", inherits = TRUE, ifnotfound = TRUE))
    list(
      L = tgd$L,
      fit = if (keep_fit) list(initializer = ini$raw, tgd = tgd,
                 k = k, rho = rho, lambda = lambda,
                 rho_rule = parameters$sgca_rho_rule[[1L]]) else NULL,
      converged = conv, completed = tgd$completed,
      status = if (tgd$completed && (!SGCA_CV_REQUIRE_CONVERGENCE || isTRUE(conv))) "ok" else "not_converged",
      iterations = init_iters + tgd$iterations,
      diagnostics = list(sgca_rho_actual = rho,
        sgca_init_iterations = init_iters, sgca_init_change = ini$convergence,
        sgca_init_converged = init_conv,
        sgca_init_stop_reason = if (isTRUE(init_conv)) "absolute_change" else "iteration_limit",
        sgca_tgd_iterations = tgd$iterations, sgca_tgd_change = tgd$absolute_change,
        sgca_tgd_converged = tgd$converged, sgca_tgd_stop_reason = tgd$stop_reason,
        sgca_stopping = SGCA_STOPPING, sgca_cv_score = SGCA_CV_SCORE,
        sgca_implementation = "paper_algorithm1_0.2.16")
    )
  }
  cross_validate_loading_grid(
    full_views, full_prep, fold_objects, rank, grid,
    prepare_context, fit_candidate, label = "SGCA", seed = seed,
    parallel_folds = parallel_folds,
    score_loss = if (SGCA_CV_SCORE == "paper") sgca_paper_validation_loss else validation_loss,
    require_convergence = SGCA_CV_REQUIRE_CONVERGENCE,
    require_completion = TRUE
  )
}
