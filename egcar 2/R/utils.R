# Input checks and scoped random-seed, worker-plan and thread utilities.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.egcar_scalar <- function(x, name, lower = 0, strict = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict) x <= lower else x < lower) ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop(name, " must be a finite ", if (integer) "integer " else "number ",
         if (strict) "> " else ">= ", lower, ".", call. = FALSE)
  }
  invisible(x)
}

.egcar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

.egcar_grid <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(if (positive) x <= 0 else x < 0))
    stop(name, " must contain finite ", if (positive) "positive" else "nonnegative",
         " numbers.", call. = FALSE)
  unique(as.numeric(x))
}

# Scope RNG and optional BLAS changes to one call, including failure paths.
.egcar_with_seed <- function(seed = NULL, code) {
  if (!is.null(seed)) .egcar_scalar(seed, "seed", 0, integer = TRUE)
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(as.integer(seed))
  force(code)
}

.egcar_with_threads <- function(threads, code) {
  if (!is.null(threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    old_omp <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
    on.exit({
      if (!is.null(old_blas)) try(RhpcBLASctl::blas_set_num_threads(old_blas), silent = TRUE)
      if (!is.null(old_omp)) try(RhpcBLASctl::omp_set_num_threads(old_omp), silent = TRUE)
    }, add = TRUE)
    try(RhpcBLASctl::blas_set_num_threads(as.integer(threads)), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(as.integer(threads)), silent = TRUE)
  }
  force(code)
}

.egcar_with_workers <- function(workers, code) {
  if (workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE))
      stop("workers > 1 requires both future and future.apply.", call. = FALSE)
    old <- future::plan()
    matching <- inherits(old, "multisession") &&
      as.integer(future::nbrOfWorkers()) == as.integer(workers)
    if (!matching) {
      on.exit(future::plan(old), add = TRUE)
      future::plan(future::multisession, workers = as.integer(workers))
    }
  }
  force(code)
}

.egcar_views <- function(views, min_rows = 2L) {
  if (!is.list(views) || is.data.frame(views) || length(views) < 2L)
    stop("views must be a list of at least two numeric matrices.", call. = FALSE)
  out <- lapply(seq_along(views), function(k) {
    X <- views[[k]]
    if (is.data.frame(X)) {
      if (!all(vapply(X, is.numeric, logical(1))))
        stop("Every data-frame column must be numeric (view ", k, ").")
      X <- as.matrix(X)
    }
    if (!is.matrix(X) || !is.numeric(X) || ncol(X) < 1L || nrow(X) < min_rows ||
        any(!is.finite(X)))
      stop("View ", k, " must be a finite numeric matrix with at least ",
           min_rows, " rows and one column.", call. = FALSE)
    storage.mode(X) <- "double"
    if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
    if (anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X)))
      stop("Column names within each view must be unique and nonempty.")
    X
  })
  n <- vapply(out, nrow, integer(1))
  if (any(n != n[[1L]])) stop("All views must have the same number of rows.")
  nm <- names(views)
  if (is.null(nm)) nm <- paste0("view", seq_along(out))
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("View names must be unique and nonempty.")
  names(out) <- nm
  out
}
