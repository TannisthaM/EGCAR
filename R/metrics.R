# Held-out generalized-Rayleigh scoring and subspace/support diagnostics.

#' Extract loadings, calculate scores and predict components
#'
#' @description Access stacked or view-wise coefficients, compute component scores with stored training means, or evaluate the common generalized Rayleigh score.
#' @param object An \code{egcar_fit} or \code{egcar_cv} object with a valid loading matrix.
#' @param newdata Complete matched views with training feature names and column order. Named views may be supplied in a different view order and are matched by name.
#' @param ridge Relative score-normalization ridge, default 1e-8.
#' @details Use \code{coef(object)} for the stacked loading matrix or \code{coef(object, by_view = TRUE)} for its view blocks. Use \code{predict(object, newdata, type = "views")} for a list of per-view component scores, or \code{type = "sum"} to sum them. Prediction does not impute an absent view. \code{egcar_score} reports the positive score; CV minimizes its negative.
#' @return A scalar validation score. The coefficient and prediction methods return matrices or lists of matrices.
#' @rdname egcar_predictions
#' @export
egcar_score <- function(object, newdata, ridge = 1e-8) {
  if (!inherits(object, "egcar_fit") && !inherits(object, "egcar_cv"))
    stop("object must be an egcar_fit or egcar_cv object.")
  .egcar_scalar(ridge, "ridge")
  # Also validates view/column order without estimating any validation means.
  .egcar_predict(object, newdata, "views")
  named <- !is.null(names(newdata))
  newdata <- .egcar_views(newdata, min_rows = 1L)
  if (named) newdata <- newdata[object$view_names]
  val <- make_validation_covariance(center_views_at(newdata, object$means))
  validation_score(object$L, val, ridge = ridge)
}

#' @rdname egcar_simulation
#' @export
egcar_distance <- function(A, B, rank = ncol(A)) {
  A <- as.matrix(A); B <- as.matrix(B)
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  if (!is.numeric(A) || !is.numeric(B) || nrow(A) != nrow(B) ||
      any(!is.finite(A)) || any(!is.finite(B))) stop("A and B must be finite numeric matrices with equal row counts.")
  sine_theta_distance(A, B, as.integer(rank))
}

validation_score <- function(L, validation, ridge = 1e-8) {
  if (is.null(L) || any(!is.finite(L))) return(-Inf)
  Q <- symmetrize(crossprod(L, validation$Sigma0 %*% L))
  scale_diag <- mean(diag(Q))
  if (!is.finite(scale_diag) || scale_diag <= 0) return(-Inf)
  ev <- eigen(Q, symmetric = TRUE)
  d <- pmax(ev$values + ridge * scale_diag, 1e-10)
  A <- crossprod(L, validation$Sigma %*% L)
  # tr(Q^(-1/2) A Q^(-1/2)) = sum_j (v_j' A v_j) / d_j.
  # No inverse matrix or two extra matrix-matrix products are constructed.
  score <- sum(colSums(ev$vectors * (A %*% ev$vectors)) / d)
  if (is.finite(score)) score else -Inf
}

validation_loss <- function(L, validation, ridge = 1e-8) {
  -validation_score(L, validation, ridge = ridge)
}

orthonormal_basis <- function(A, tol = 1e-10) {
  ss <- svd(A, nu = min(dim(A)), nv = 0)
  if (length(ss$d) == 0L || ss$d[[1L]] <= 0) return(NULL)
  keep <- which(ss$d > tol * ss$d[[1L]])
  if (length(keep) == 0L) return(NULL)
  ss$u[, keep, drop = FALSE]
}

sine_theta_distance <- function(A, B, rank) {
  QA <- orthonormal_basis(A)
  QB <- orthonormal_basis(B)
  if (is.null(QA) || is.null(QB) || ncol(QA) < rank || ncol(QB) < rank) {
    return(NA_real_)
  }
  QA <- QA[, seq_len(rank), drop = FALSE]
  QB <- QB[, seq_len(rank), drop = FALSE]
  cc <- svd(crossprod(QA, QB), nu = 0, nv = 0)$d
  cc <- pmin(1, pmax(0, cc))
  sqrt(sum(pmax(0, 1 - cc^2)))
}

support_metrics <- function(C_full, true_active, threshold = 1e-6) {
  selected <- which(row_l2(C_full) > threshold)
  tp <- length(intersect(selected, true_active))
  fp <- length(setdiff(selected, true_active))
  fn <- length(setdiff(true_active, selected))
  precision <- if (length(selected) == 0L) 0 else tp / length(selected)
  recall <- if (length(true_active) == 0L) 1 else tp / length(true_active)
  fdp <- if (length(selected) == 0L) 0 else fp / length(selected)
  c(precision = precision, recall = recall, fdp = fdp, tp = tp, fp = fp, fn = fn)
}
