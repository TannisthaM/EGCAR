# Accelerated EGCAR: public fit, rate, regularization-path and CV entry points.

#' Fit accelerated entrywise-only or group-only EGCAR
#'
#' @description Fit the pairwise regression loss with one of two separate sparsity penalties and obtain localized generalized-correlation loadings.
#' @param x Raw views, an \code{egcar_prepared} object, or an \code{egcar_cv_data} object (using its full data).
#' @param rank Number of leading positive loading directions requested.
#' @param penalty \code{"l11"} for entrywise sparsity or \code{"l21"} for the view-wise concatenated row-group penalty. They are not combined.
#' @param lambda A direct nonnegative penalty coefficient, or the coefficient vector for a path.
#' @param control An \code{egcar_control} object or a named subset of its arguments.
#' @param init An optional same-family fitted object or dimension-compatible ADMM-state list. Used as a warm start.
#' @param multiplier Nonnegative multiplier applied to the selected rate scale.
#' @details For each unordered view pair the smooth loss is
#' \deqn{\tfrac12\mathrm{tr}(C_{kl}^{T}S_{kk}C_{kl}S_{ll})-\langle S_{kl},C_{kl}\rangle.}
#' The entrywise version adds \eqn{\lambda\sum_{k<l}\|C_{kl}\|_{1,1}}.
#' The group version adds \eqn{\lambda\sum_k\|M_k(C)\|_{2,1}}, where each row
#' of \eqn{M_k(C)} concatenates all coefficients incident to one variable in view k.
#' The latter uses only C, G and V; it contains no entrywise proximal update.
#' The reduced-basis linear solve retains the full null-space contribution.
#' Loadings use the largest algebraic eigenvalues of the localized normalized operator.
#' If fewer than the requested number of positive directions exist, the operator fit
#' is returned with status \code{"invalid_loading"}; no lower rank is silently substituted.
#' The rate scales are \eqn{\sqrt{\log(p)/n}} and
#' \eqn{\sqrt{(d_{max}+\log(p))/n}}, respectively, where
#' \eqn{d_{max}=\max_k(p-p_k)}. The group scale is inherited as a numerical
#' benchmark normalization, not asserted to be an optimal theoretical rate.
#' Paths visit coefficients from strongest to weakest, return them in input order,
#' and reuse primal/scaled-dual states with the configured initial augmentation
#' parameter, matching the supplied script's convention.
#' @return \code{egcar_fit} and \code{egcar_rate} return \code{egcar_fit} objects containing stacked loadings \code{L}, view loadings, edge blocks \code{C}, the raw solver state, residuals, convergence status, selected rows, training means, controls and separated timings. A path returns the coefficient vector and corresponding fitted objects.
#' @rdname egcar_fit
#' @examples
#' sim <- egcar_simulate(n = 30, p_list = c(3, 4, 5), active_per_view = 2)
#' fit <- egcar_fit(sim$views, penalty = "l21", lambda = 0.01)
#' print(fit)
#' @export
egcar_fit <- function(x, rank = 1L, penalty = c("l11", "l21"), lambda = 0.01,
                      control = egcar_control(), init = NULL) {
  cl <- match.call()
  penalty <- match.arg(penalty)
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  .egcar_scalar(lambda, "lambda")
  control <- .egcar_as_control(control)
  .egcar_with_seed(NULL, .egcar_with_threads(control$blas_threads, {
    prepared <- egcar_prepare(x)
    if (rank > prepared$p) stop("rank exceeds the total number of variables.")
    e <- .egcar_engine(control)
    prep <- prepared$prep
    prep$loading_factor_cache <- new.env(parent = emptyenv())
    t0 <- proc.time()[[3L]]
    raw <- .egcar_solve_one(prep, penalty, lambda, control, e, init)
    fit_time <- proc.time()[[3L]] - t0
    t0 <- proc.time()[[3L]]
    loading <- .egcar_loading(prep, raw, as.integer(rank), control, e)
    loading_time <- proc.time()[[3L]] - t0
    .egcar_fit_object(raw, loading, prepared, as.integer(rank), penalty, lambda,
                      control, fit_time, loading_time, cl)
  }))
}

#' @rdname egcar_fit
#' @export
egcar_rate <- function(x, rank = 1L, penalty = c("l11", "l21"), multiplier = 1,
                       control = egcar_control()) {
  penalty <- match.arg(penalty)
  .egcar_scalar(multiplier, "multiplier")
  prepared <- egcar_prepare(x)
  if (!is.finite(prepared$n)) stop("Rate scaling requires a finite sample size.")
  lambda <- multiplier * if (penalty == "l11") sqrt(log(prepared$p) / prepared$n) else
    sqrt((max(prepared$p - prepared$p_list) + log(prepared$p)) / prepared$n)
  out <- egcar_fit(prepared, rank, penalty, lambda, control)
  out$rate_multiplier <- multiplier
  out$rate_rule <- if (penalty == "l11") "sqrt(log(p)/n)" else "sqrt((d_max+log(p))/n)"
  out$call <- match.call()
  out
}

#' @rdname egcar_fit
#' @export
egcar_path <- function(x, rank = 1L, penalty = c("l11", "l21"),
                       lambda = 10^seq(-5, 4), control = egcar_control()) {
  penalty <- match.arg(penalty)
  lambda <- .egcar_grid(lambda, "lambda")
  prepared <- egcar_prepare(x)
  fits <- vector("list", length(lambda))
  previous <- NULL
  for (i in order(lambda, decreasing = TRUE)) {
    fits[[i]] <- egcar_fit(prepared, rank, penalty, lambda[[i]], control, init = previous)
    previous <- fits[[i]]
  }
  names(fits) <- format(lambda, digits = 17)
  structure(list(lambda = lambda, fits = fits, penalty = penalty, rank = rank,
                 traversal = order(lambda, decreasing = TRUE)), class = "egcar_path")
}

#' Cross-validate either accelerated EGCAR penalty
#'
#' @description Evaluate a full one-dimensional direct-coefficient grid on common folds, select the highest average held-out generalized Rayleigh score, and refit on all observations.
#' @param x Raw views or a reusable \code{egcar_cv_data} object.
#' @param rank Requested positive loading rank.
#' @param penalty Separate entrywise \code{"l11"} or row-group \code{"l21"} estimator.
#' @param lambda Direct coefficient grid, strictly positive; zero is excluded. No fold-size scaling is applied.
#' @param fold_id Optional supplied fold labels, or NULL to use labels already stored in x.
#' @param nfolds Number of folds when creating them from raw views.
#' @param seed Nonnegative integer seed; the caller's RNG state is restored.
#' @param workers Fold-level workers, capped at the number of folds. Values above one require future and future.apply. The previous future plan is restored.
#' @param control An \code{egcar_control} object or named argument list.
#' @details Candidates are visited in descending coefficient order within each fold
#' and use same-family warm starts. Finite but nonconverged fits remain eligible,
#' as in the supplied benchmark; convergence flags and iteration counts are recorded.
#' Every fold must have a finite score for a candidate to be eligible. Exact score
#' ties favor stronger regularization. Zero is not an admissible CV candidate for either EGCAR family.
#'
#' The common loss is the negative generalized Rayleigh score, with the supplied
#' relative ridge and eigenvalue floor. It is not reconstruction MSE. Validation
#' covariances use training-mean centering only. If every candidate fails, the
#' result has status \code{"no_valid_cv"}, no selected penalty and no refit. This
#' explicit failure replaces the earlier script's fallback to an arbitrary
#' coefficient of one. Diagnostics are in \code{cv_fold_table}.
#'
#' Preparation is reported separately and can be amortized over methods. Tuning
#' time includes fold fits and scoring, and may include worker startup. EGCAR's
#' final \code{fit_time} measures ADMM, \code{loading_time} its final loading
#' extraction, and \code{total_time} their sum plus tuning. The time of loading
#' or compiling the installed package is not included.
#' @return An \code{egcar_cv} object containing \code{fit}, \code{L}, \code{best}, \code{lambda}, candidate and fold tables, supplied fold labels, status, controls and timings.
#' @rdname egcar_cv
#' @examples
#' sim <- egcar_simulate(n = 24, p_list = c(3, 4, 5), active_per_view = 2)
#' shared <- egcar_cv_data(sim$views, nfolds = 2)
#' cv <- egcar_cv(shared, penalty = "l11", lambda = c(0.01, 0.1))
#' print(cv)
#' @export
egcar_cv <- function(x, rank = 1L, penalty = c("l11", "l21"),
                     lambda = 10^seq(-5, 4), fold_id = NULL,
                     nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control()) {
  cl <- match.call()
  penalty <- match.arg(penalty)
  lambda <- egcar_positive_cv_grid(lambda, "lambda")
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  .egcar_scalar(workers, "workers", 1, integer = TRUE)
  control <- .egcar_as_control(control)
  data <- .egcar_cv_input(x, fold_id, nfolds, seed)
  if (rank > data$p) stop("rank exceeds the number of variables.")
  workers <- as.integer(min(workers, data$nfolds))
  e <- .egcar_engine(control, workers = workers)
  order_path <- order(lambda, decreasing = TRUE)
  timing_start <- proc.time()[[3L]]
  .egcar_with_seed(seed, .egcar_with_threads(control$blas_threads,
    .egcar_with_workers(workers, {
      rows <- e$parallel_map_candidates(data$folds, function(fo) {
        f <- fo$fold
        fo$prep$loading_factor_cache <- new.env(parent = emptyenv())
        previous <- NULL
        result <- vector("list", length(lambda))
        for (j in order_path) {
          start <- proc.time()[[3L]]
          notes <- character()
          one <- tryCatch(withCallingHandlers({
            raw <- .egcar_solve_one(fo$prep, penalty, lambda[[j]], control, e,
                                   init = previous, cv = TRUE, history = FALSE)
            loading <- .egcar_loading(fo$prep, raw, as.integer(rank), control, e, FALSE)
            score <- if (isTRUE(loading$valid)) validation_score(loading$L, fo$validation) else -Inf
            list(score = score, converged = raw$converged, iterations = raw$iterations,
                 state = raw[if (penalty == "l11") c("C", "Z", "H") else c("C", "G", "V")],
                 error = if (!isTRUE(loading$valid)) loading$reason else
                   if (!is.finite(score)) "Nonfinite common validation score." else NA_character_)
          }, warning = function(w) {
            notes <<- unique(c(notes, conditionMessage(w)))
            invokeRestart("muffleWarning")
          }), error = function(err) {
            list(score = -Inf, converged = FALSE, iterations = NA_integer_,
                 state = NULL, error = conditionMessage(err))
          })
          if (!is.null(one$state)) previous <- one$state
          result[[j]] <- data.frame(candidate = j, lambda = lambda[[j]], fold = f,
            n_train = fo$prep$n, score = one$score, loss = -one$score,
            converged = one$converged, iterations = one$iterations,
            elapsed = proc.time()[[3L]] - start, error_message = one$error,
            warning_message = if (length(notes)) paste(notes, collapse = " | ") else NA_character_)
        }
        do.call(rbind, result)
      })
      fold_table <- do.call(rbind, rows)
      grid <- data.frame(lambda = lambda, candidate = seq_along(lambda))
      table <- summarize_loading_cv(grid, fold_table, data$nfolds)
      valid <- which(is.finite(table$mean_loss))
      best_index <- if (length(valid)) {
        valid[order(table$mean_loss[valid], -table$lambda[valid])[[1L]]]
      } else NA_integer_
      table$selected <- !is.na(best_index) & seq_len(nrow(table)) == best_index
      tuning_time <- proc.time()[[3L]] - timing_start
      out <- list(method = paste0("EGCAR-", toupper(penalty)), penalty = penalty,
        rank = as.integer(rank), cv_table = table, cv_fold_table = fold_table,
        fold_id = data$fold_id, fold_labels = data$fold_labels, workers = workers,
        mean_loss = if (length(valid)) table$mean_loss[[best_index]] else Inf,
        tuning_time = tuning_time, preparation_time = data$preparation_time,
        means = data$full$means, view_names = data$full$view_names,
        feature_names = data$full$feature_names, p_list = data$p_list,
        p = data$p, n = data$n, call = cl, control = control)
      if (is.na(best_index)) {
        out$fit <- NULL; out$L <- NULL; out$best <- NULL; out$lambda <- NA_real_
        out$fit_time <- 0; out$loading_time <- 0; out$total_time <- tuning_time
        out$status <- "no_valid_cv"; out$converged <- FALSE
        out$iterations <- NA_integer_
        out$error <- "No candidate had a finite requested-rank score on every fold. Inspect cv_fold_table."
      } else {
        out$best <- table[best_index, , drop = FALSE]
        out$lambda <- table$lambda[[best_index]]
        refit_start <- proc.time()[[3L]]
        refit <- tryCatch(egcar_fit(data$full, rank, penalty, out$lambda, control),
                          error = function(err) list(error = conditionMessage(err)))
        if (!inherits(refit, "egcar_fit")) {
          out$fit <- NULL; out$L <- NULL
          out$fit_time <- proc.time()[[3L]] - refit_start; out$loading_time <- 0
          out$status <- "refit_error"; out$converged <- FALSE
          out$iterations <- NA_integer_; out$error <- refit$error
        } else {
          out$fit <- refit; out$L <- refit$L; out$fit_time <- refit$fit_time
          out$loading_time <- refit$loading_time; out$status <- refit$status
          out$converged <- refit$converged; out$iterations <- refit$iterations
          out$error <- refit$error
        }
        out$total_time <- tuning_time + out$fit_time + out$loading_time
      }
      class(out) <- "egcar_cv"
      out
    })))
}
