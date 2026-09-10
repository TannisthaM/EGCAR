library(egcar)
# Long optional tests are opt-in; nothing is silently installed by a test.
if (identical(Sys.getenv("EGCAR_TEST_OPTIONAL", "false"), "true")) {
  sim <- egcar_simulate(n = 40, p_list = c(4, 5, 6), rank = 1,
                        active_per_view = 3, seed = 85)
  shared <- egcar_cv_data(sim$views, nfolds = 2, seed = 86)
  bc <- benchmark_control(sgca_init_max_iter = 60, sgca_tgd_max_iter = 100,
                           rgcca_max_iter = 200, multicca_niter = 10)
  ct <- egcar_control(max_iter_cv = 150, max_iter = 300)
  check <- function(fit) {
    stopifnot(identical(fit$fold_id, shared$fold_id), nrow(fit$cv_fold_table) >= 2)
    stopifnot(all(fit$cv_fold_table$n_train == 20))
    if (!is.null(fit$L)) stopifnot(identical(dim(fit$L), c(15L, 1L)))
  }
  if (requireNamespace("RGCCA", quietly = TRUE)) {
    a <- rgcca_cv(shared, tau_grid = c(0.1, 1), benchmarks = bc); check(a)
    b <- sgcca_cv(shared, sparsity_grid = c(0.7, 1), benchmarks = bc); check(b)
    stopifnot(!is.null(a$best), !is.null(b$best))
  }
  if (requireNamespace("PMA", quietly = TRUE)) {
    a <- multicca_cv(shared, penalty_grid = c(1.2, 2), benchmarks = bc); check(a)
    stopifnot(!is.null(a$best))
    # Compare matrix-only and public PMA paths on identical folds and starts.
    bc$multicca_backend <- "PMA"
    b <- multicca_cv(shared, penalty_grid = c(1.2, 2), benchmarks = bc); check(b)
    stopifnot(isTRUE(all.equal(a$cv_table$mean_loss, b$cv_table$mean_loss, tolerance = 1e-5)))
  }
  {
    a <- sgca_cv(shared, k_grid = c(6, 10), rho_grid = 0.01,
                  lambda_grid = 0.1, benchmarks = bc); check(a)
  }
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)) {
    old <- future::plan()
    a <- egcar_cv(shared, 1, "l21", lambda = c(0.02, 0.01), workers = 1, control = ct)
    b <- egcar_cv(shared, 1, "l21", lambda = c(0.02, 0.01), workers = 2, control = ct)
    stopifnot(isTRUE(all.equal(a$cv_table$mean_loss, b$cv_table$mean_loss, tolerance = 1e-8)))
    stopifnot(identical(class(old), class(future::plan())))
  }
}
