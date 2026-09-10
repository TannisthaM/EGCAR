# SGCCA comparison: common-loss CV over view-wise sparsity bounds.

#' @rdname comparison_cv
#' @export
sgcca_cv <- function(x, rank = 1L, sparsity_grid = NULL, fold_id = NULL,
                     nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control(), benchmarks = benchmark_control()) {
  .egcar_external_cv(x, rank, "SGCCA", fold_id, nfolds, seed, workers, control,
                     benchmarks, list(sparsity_grid = sparsity_grid), match.call())
}

sgcca_common_cv <- function(
    full_views, full_prep, fold_objects, rank,
    sparsity_grid = SGCCA_SPARSITY_GRID, seed = 1L, parallel_folds = PARALLEL_CV) {
  rgcca_family_common_cv(
    full_views, full_prep, fold_objects, rank, sparsity_grid,
    method = "sgcca", seed = seed, parallel_folds = parallel_folds
  )
}
