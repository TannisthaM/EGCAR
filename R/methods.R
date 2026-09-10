# Print, coefficient and prediction methods for package result classes.

#' @rdname egcar_prepare
#' @export
#' @param x A fitted or prepared object of the corresponding class.
#' @param ... Further arguments reserved for S3 compatibility.
print.egcar_prepared <- function(x, ...) {
  cat("Prepared EGCAR data:", x$n, "observations; view dimensions",
      paste(x$p_list, collapse = ", "), "\n")
  invisible(x)
}

#' @rdname egcar_prepare
#' @export
#' @param x A fitted or prepared object of the corresponding class.
#' @param ... Further arguments reserved for S3 compatibility.
print.egcar_cv_data <- function(x, ...) {
  cat("Shared EGCAR CV data:", x$n, "observations;", x$nfolds, "folds; dimensions",
      paste(x$p_list, collapse = ", "), "\n")
  cat("Training-mean centering only; covariance divisor is the split sample size.\n")
  invisible(x)
}

#' @rdname egcar_fit
#' @export
#' @param x A fitted or prepared object of the corresponding class.
#' @param ... Further arguments reserved for S3 compatibility.
print.egcar_fit <- function(x, ...) {
  cat("Accelerated EGCAR:", x$penalty, "only; lambda =", format(x$lambda),
      "; rank =", x$rank, "\n")
  cat("Status:", x$status, "; ADMM iterations:", x$iterations,
      "; selected loading rows:", length(x$selected), "\n")
  if (!is.null(x$error)) cat("Diagnostic:", x$error, "\n")
  invisible(x)
}

#' @rdname egcar_cv
#' @export
#' @param x A fitted or prepared object of the corresponding class.
#' @param ... Further arguments reserved for S3 compatibility.
print.egcar_cv <- function(x, ...) {
  cat(x$method, "common-loss CV; rank =", x$rank, "; status:", x$status, "\n")
  if (!is.null(x$best)) {
    cols <- intersect(c("lambda", "sgca_k", "sgca_rho", "sgca_lambda", "rgcca_tau",
                        "sgcca_sparsity", "multicca_l1_bound", "mean_loss"), names(x$best))
    print(x$best[, cols, drop = FALSE], row.names = FALSE)
  } else cat(x$error, "\n")
  cat("Tuning:", format(x$tuning_time, digits = 4), "s; final refit:",
      format(x$fit_time + x$loading_time, digits = 4), "s\n")
  invisible(x)
}

#' @rdname egcar_predictions
#' @export
#' @param by_view Return one loading matrix per view rather than the stacked matrix.
#' @param ... Further arguments reserved for S3 compatibility.
coef.egcar_fit <- function(object, by_view = FALSE, ...) {
  .egcar_flag(by_view, "by_view")
  if (is.null(object$L)) stop("This fit has no valid requested-rank loading: ", object$status)
  if (by_view) object$loadings else object$L
}

#' @rdname egcar_predictions
#' @export
#' @param by_view Return one loading matrix per view rather than the stacked matrix.
#' @param ... Further arguments reserved for S3 compatibility.
coef.egcar_cv <- function(object, by_view = FALSE, ...) {
  .egcar_flag(by_view, "by_view")
  if (is.null(object$L)) stop("This CV result has no valid loading: ", object$status)
  if (!by_view) return(object$L)
  setNames(lapply(make_block_indices(object$p_list), function(ii) object$L[ii, , drop = FALSE]),
           object$view_names)
}

.egcar_predict <- function(object, newdata, type) {
  if (is.null(object$L)) stop("This object has no valid loading matrix.")
  type <- match.arg(type, c("views", "sum"))
  named <- !is.null(names(newdata))
  newdata <- .egcar_views(newdata, min_rows = 1L)
  if (length(newdata) != length(object$p_list)) stop("newdata has a different number of views.")
  if (named) {
    if (!setequal(names(newdata), object$view_names)) stop("newdata view names do not match the fit.")
    newdata <- newdata[object$view_names]
  }
  idx <- make_block_indices(object$p_list)
  scores <- lapply(seq_along(newdata), function(k) {
    X <- newdata[[k]]
    if (ncol(X) != object$p_list[[k]]) stop("newdata dimensions do not match in view ", k)
    # No implicit column reordering: all columns must retain training order.
    if (!identical(colnames(X), object$feature_names[[k]]))
      stop("newdata columns must have the same names and order as training data in view ", k)
    sweep(X, 2L, object$means[[k]], "-") %*% object$L[idx[[k]], , drop = FALSE]
  })
  names(scores) <- object$view_names
  if (type == "sum") Reduce(`+`, scores) else scores
}

#' @rdname egcar_predictions
#' @export
#' @param type Return per-view component scores ("views") or their sum ("sum").
#' @param ... Further arguments reserved for S3 compatibility.
predict.egcar_fit <- function(object, newdata, type = c("views", "sum"), ...) {
  .egcar_predict(object, newdata, match.arg(type))
}

#' @rdname egcar_predictions
#' @export
#' @param type Return per-view component scores ("views") or their sum ("sum").
#' @param ... Further arguments reserved for S3 compatibility.
predict.egcar_cv <- function(object, newdata, type = c("views", "sum"), ...) {
  .egcar_predict(object, newdata, match.arg(type))
}
