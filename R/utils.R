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

.egcar_optional_limit <- function(x, name, integer = FALSE) {
  if (is.numeric(x) && length(x) == 1L && !is.na(x) && identical(as.double(x), Inf))
    return(invisible(x))
  .egcar_scalar(x, name, if (integer) 1 else 0, strict = !integer, integer = integer)
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

.egcar_worker_plan_is_multicore <- function() {
  # future::multicore forks the current R process, so each worker shares the
  # parent's already-computed, read-only data (prepared covariance blocks,
  # cached eigenbases, etc.) via copy-on-write -- nothing is duplicated in
  # physical memory unless a worker actually writes to it. future::multisession
  # instead spawns fully independent R sessions with no shared memory at all,
  # so every worker must hold its own complete copy of whatever it touches.
  # At large p (many GB of working set per worker) this is the difference
  # between ~1x and ~workers-x peak memory. multicore requires OS-level
  # forking (Unix-like only, and not always available even there, e.g. some
  # restricted interactive front-ends), so this checks support rather than
  # just assuming it, and multisession remains the safe, portable fallback
  # everywhere it isn't. Set EGCAR_MULTICORE=0 to force multisession even
  # where multicore would otherwise be used.
  if (identical(Sys.getenv("EGCAR_MULTICORE", "1"), "0")) return(FALSE)
  .Platform$OS.type == "unix" &&
    requireNamespace("parallelly", quietly = TRUE) &&
    isTRUE(tryCatch(parallelly::supportsMulticore(), error = function(e) FALSE))
}

.egcar_with_workers <- function(workers, code) {
  if (workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE))
      stop("workers > 1 requires both future and future.apply.", call. = FALSE)
    use_multicore <- .egcar_worker_plan_is_multicore()
    plan_class <- if (use_multicore) "multicore" else "multisession"
    old <- future::plan()
    matching <- inherits(old, plan_class) &&
      as.integer(future::nbrOfWorkers()) == as.integer(workers)
    if (!matching) {
      on.exit(future::plan(old), add = TRUE)
      # future::plan() substitute()s its strategy argument by default, so it
      # must be a plain value/symbol, never an inline if/else expression --
      # passing the conditional directly here (as the 0.2.8 patch originally
      # did) made plan() capture the literal call `if(use_multicore,
      # future::multicore, future::multisession)` instead of picking one of
      # the two functions first. Its internal tweak() then tried to treat
      # the primitive `if` itself as the strategy and inspect its
      # environment -- primitives have none (environment(`if`) is NULL) --
      # crashing with "Error in ls(envir = env, ...) : invalid 'envir'
      # argument" inside tweak.function, unconditionally, every time,
      # regardless of e/dispatcher/globals (confirmed via a full traceback:
      # the crash is in plan() setup, before any CV-fold work is reached).
      # Evaluating the choice into a plain variable first, and passing that
      # to plan(), gives it a simple symbol to substitute() instead.
      strategy_fn <- if (use_multicore) future::multicore else future::multisession
      future::plan(strategy_fn, workers = as.integer(workers))
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

# Normalize solver state at the R boundary. This additionally supports legacy
# nested Rcpp output in which an arma::mat lost only its dim attribute. The
# dimensions are not guessed: each edge has a unique known p_k x p_l shape.
.egcar_normalize_solver_state <- function(state, context, group) {
  fields <- if (group) {
    if (!is.null(state[["a"]])) c("C", "Hk", "Hl") else c("C", "Gk", "Gl", "Vk", "Vl")
  } else c("C", "Z", "H")
  if (!is.list(state)) return(state)
  count <- length(context$edge_k)
  for (nm in intersect(fields, names(state))) {
    blocks <- state[[nm]]
    if (!is.list(blocks) || length(blocks) != count) next
    state[[nm]] <- lapply(seq_len(count), function(j) {
      k <- context$edge_k[[j]]; l <- context$edge_l[[j]]
      nr <- as.integer(context$p_list[[k]]); nc <- as.integer(context$p_list[[l]])
      A <- blocks[[j]]
      if (is.matrix(A) && is.numeric(A) && identical(dim(A), c(nr, nc))) return(A)
      if (is.numeric(A) && is.null(dim(A)) && length(A) == nr * nc)
        return(matrix(as.numeric(A), nrow = nr, ncol = nc))
      A
    })
  }
  state
}

# Validate the native/optimized solver boundary before loading extraction.
.egcar_check_solver_state <- function(state, context, group, backend) {
  fields <- if (group) {
    if (!is.null(state[["a"]])) c("C", "Hk", "Hl") else c("C", "Gk", "Gl", "Vk", "Vl")
  } else c("C", "Z", "H")
  hint <- if (identical(backend, "cpp")) paste0(
    " Reinstall the patched egcar source package and restart R, including CV workers.") else ""
  fail <- function(message) stop("EGCAR ", backend, " matrix interface: ", message,
                                  hint, call. = FALSE)
  if (!is.list(state) || anyDuplicated(names(state)) ||
      !all(fields %in% names(state)))
    fail("the solver did not return the required named state lists.")
  count <- length(context$edge_k)
  if (group && !is.null(state[["a"]])) {
    if (!is.list(state$a) || length(state$a) != length(context$p_list))
      fail("invalid group multiplier list a.")
    for (k in seq_along(context$p_list))
      if (!is.numeric(state$a[[k]]) || !is.null(dim(state$a[[k]])) ||
          length(state$a[[k]]) != context$p_list[[k]] ||
          any(!is.finite(state$a[[k]])) || any(state$a[[k]] < 0 | state$a[[k]] > 1))
        fail("invalid group row multipliers.")
  }
  for (nm in fields) {
    blocks <- state[[nm]]
    if (!is.list(blocks) || length(blocks) != count)
      fail(paste0("state$", nm, " must contain ", count, " edge matrices."))
    for (j in seq_len(count)) {
      expected <- as.integer(c(context$p_list[[context$edge_k[[j]]]],
                               context$p_list[[context$edge_l[[j]]]]))
      A <- blocks[[j]]
      if (!is.matrix(A) || !is.numeric(A) || !identical(dim(A), expected)) {
        observed <- if (is.null(dim(A))) paste0(
          "a dimensionless ", typeof(A), " object of length ", length(A)) else paste0(
          "a ", typeof(A), " array with dimensions ", paste(dim(A), collapse = " x "))
        fail(paste0("state$", nm, "[[", j, "]] (edge ",
          context$edge_k[[j]], "_", context$edge_l[[j]], ") must be a numeric ",
          paste(expected, collapse = " x "), " matrix; got ", observed, "."))
      }
    }
  }
  invisible(TRUE)
}
