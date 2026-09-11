# Data preparation and reusable, training-centered cross-validation folds.

# Cached statistics are independent of the regularization coefficient.
#' Prepare views and reusable common cross-validation folds
#'
#' @description Compute cached covariance blocks, spectral reductions, training means and, when requested, common folds for all methods.
#' @param views A list of at least two finite numeric matrices with matching sample rows. Numeric data frames are accepted. Variable and block scaling are not performed. Missing observations and missing views are not supported.
#' @param n Number of observations.
#' @param nfolds Number of nonempty folds.
#' @param seed Nonnegative integer seed. The caller's RNG state is restored.
#' @param fold_id One nonmissing fold label per row. Supplied labels are never regenerated.
#' @details Every covariance uses divisor n for that split. Validation observations are centered at training means, not their own means. All training splits must have at least two observations. A prepared object stores covariance information and means but not the full raw data. A CV-data object retains raw data and full and fold preparations, so it can be reused for both EGCAR penalties and all comparison methods. Stored means and raw data may be sensitive; save these objects accordingly.
#' @return \code{egcar_prepare} returns an \code{egcar_prepared} object. \code{egcar_folds} returns integer fold labels. \code{egcar_cv_data} returns an \code{egcar_cv_data} object.
#' @rdname egcar_prepare
#' @examples
#' sim <- egcar_simulate(n = 30, p_list = c(4, 5, 6), active_per_view = 3)
#' shared <- egcar_cv_data(sim$views, nfolds = 3)
#' print(shared)
#' @export
egcar_prepare <- function(views) {
  if (inherits(views, "egcar_prepared")) return(views)
  if (inherits(views, "egcar_cv_data")) return(views$full)
  t0 <- proc.time()[[3L]]
  views <- .egcar_views(views)
  centered <- center_views(views)
  e <- .egcar_engine()
  prep <- e$prepare_problem(centered$views)
  structure(list(prep = prep, means = centered$means,
                 p_list = prep$p_list, p = prep$p, n = prep$n,
                 view_names = names(views), feature_names = lapply(views, colnames),
                 preparation_time = proc.time()[[3L]] - t0), class = "egcar_prepared")
}

#' @rdname egcar_prepare
#' @export
egcar_folds <- function(n, nfolds = 5L, seed = 1L) {
  .egcar_scalar(n, "n", 3, integer = TRUE)
  .egcar_scalar(nfolds, "nfolds", 2, integer = TRUE)
  if (nfolds > n) stop("nfolds cannot exceed n.")
  .egcar_with_seed(seed, sample(rep(seq_len(nfolds), length.out = n)))
}

#' @rdname egcar_prepare
#' @export
egcar_cv_data <- function(views, fold_id = NULL, nfolds = 5L, seed = 1L) {
  if (inherits(views, "egcar_cv_data")) {
    if (!is.null(fold_id) && !identical(as.character(fold_id), as.character(views$fold_id)))
      stop("fold_id conflicts with the supplied egcar_cv_data object.")
    return(views)
  }
  t0 <- proc.time()[[3L]]
  views <- .egcar_views(views, min_rows = 3L)
  n <- nrow(views[[1L]])
  if (is.null(fold_id)) fold_id <- egcar_folds(n, nfolds, seed)
  if ((!is.numeric(fold_id) && !is.character(fold_id) && !is.factor(fold_id)) ||
      length(fold_id) != n || anyNA(fold_id))
    stop("fold_id must give one nonmissing fold label per row.")
  labs <- sort(unique(fold_id))
  if (length(labs) < 2L) stop("At least two folds are required.")
  counts <- tabulate(match(fold_id, labs), nbins = length(labs))
  if (any(n - counts < 2L)) stop("Every training split must contain at least two rows.")
  full <- egcar_prepare(views)
  e <- .egcar_engine()
  objects <- e$make_fold_objects(views, fold_id)
  structure(list(views = views, full = full, folds = objects, fold_id = fold_id,
                 fold_labels = labs, nfolds = length(labs), n = n,
                 p_list = full$p_list, p = full$p,
                 preparation_time = proc.time()[[3L]] - t0), class = "egcar_cv_data")
}

.egcar_cv_input <- function(x, fold_id, nfolds, seed) {
  if (inherits(x, "egcar_prepared"))
    stop("Cross-validation needs raw views or egcar_cv_data(), not covariance-only egcar_prepare().")
  egcar_cv_data(x, fold_id = fold_id, nfolds = nfolds, seed = seed)
}

make_folds <- function(n, K = 5L, seed = 1L) {
  set.seed(seed)
  sample(rep(seq_len(K), length.out = n))
}

make_fold_objects <- function(views, fold_id) {
  fold_labels <- sort(unique(fold_id))
  lapply(fold_labels, function(f) {
    train_idx <- which(fold_id != f)
    val_idx <- which(fold_id == f)
    train_raw <- lapply(views, function(X) X[train_idx, , drop = FALSE])
    val_raw <- lapply(views, function(X) X[val_idx, , drop = FALSE])
    centered <- center_views(train_raw)
    train_views <- centered$views
    val_views <- center_views_at(val_raw, centered$means)
    list(
      # Keep only fields actually consumed by CV. Raw indices and means are
      # dropped after centering to reduce master/worker payloads.
      fold = f,
      train_views = train_views,
      prep = prepare_problem(train_views),
      validation = make_validation_covariance(val_views)
    )
  })
}
