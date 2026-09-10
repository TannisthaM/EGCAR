# Explicit single-penalty EGCAR entry points. No combined-penalty aliases.
#' @rdname EGCAR_methods
#' @export
EGCAR_L11 <- function(x, rank = 1L, lambda = 0.01, control = egcar_control(), init = NULL) {
  egcar_fit(x, rank = rank, penalty = "l11", lambda = lambda, control = control, init = init)
}
#' @rdname EGCAR_methods
#' @export
EGCAR_L21 <- function(x, rank = 1L, lambda = 0.01, control = egcar_control(), init = NULL) {
  egcar_fit(x, rank = rank, penalty = "l21", lambda = lambda, control = control, init = init)
}
#' @rdname EGCAR_methods
#' @export
EGCAR_L11_CV <- function(x, rank = 1L, lambda = 10^seq(-5, 4), fold_id = NULL,
    nfolds = 5L, seed = 1L, workers = 1L, control = egcar_control()) {
  egcar_cv(x, rank = rank, penalty = "l11", lambda = lambda, fold_id = fold_id,
    nfolds = nfolds, seed = seed, workers = workers, control = control)
}
#' @rdname EGCAR_methods
#' @export
EGCAR_L21_CV <- function(x, rank = 1L, lambda = 10^seq(-5, 4), fold_id = NULL,
    nfolds = 5L, seed = 1L, workers = 1L, control = egcar_control()) {
  egcar_cv(x, rank = rank, penalty = "l21", lambda = lambda, fold_id = fold_id,
    nfolds = nfolds, seed = seed, workers = workers, control = control)
}
#' @rdname EGCAR_methods
#' @export
EGCAR_L11_Rate <- function(x, rank = 1L, multiplier = 1, control = egcar_control()) {
  egcar_rate(x, rank = rank, penalty = "l11", multiplier = multiplier, control = control)
}
#' @rdname EGCAR_methods
#' @export
EGCAR_L21_Rate <- function(x, rank = 1L, multiplier = 1, control = egcar_control()) {
  egcar_rate(x, rank = rank, penalty = "l21", multiplier = multiplier, control = control)
}
