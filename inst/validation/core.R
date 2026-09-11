library(egcar)
close <- function(x, y, tolerance = 1e-6) {
  stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance, check.attributes = FALSE)))
}
assert_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

# Controls and invalid inputs should fail before entering native code.
assert_error(egcar_control(mu = 0))
assert_error(egcar_control(check_every = 0))
assert_error(egcar_control(max_iter = 1.5))
assert_error(egcar_prepare(list(matrix(1, 3, 2))))
assert_error(egcar_prepare(list(matrix(NA_real_, 3, 2), matrix(1, 3, 2))))
assert_error(egcar_cv_data(list(matrix(1, 3, 2), matrix(1, 3, 2)), fold_id = c(1, 1, 2)))

set.seed(991)
before <- .Random.seed
sim <- egcar_simulate(n = 35, p_list = c(4, 5, 6), rank = 2,
                      active_per_view = 3, seed = 15)
stopifnot(identical(before, .Random.seed))
views <- sim$views
prepared <- egcar_prepare(views)
stopifnot(inherits(prepared, "egcar_prepared"), prepared$n == 35)

# Native, optimized R, and dense-reference iterations agree, including warm starts.
for (wide in c(FALSE, TRUE)) {
  X <- if (wide) list(matrix(rnorm(15 * 20), 15, 20),
                       matrix(rnorm(15 * 8), 15, 8), matrix(rnorm(15 * 5), 15, 5)) else views
  X[[1]][, 1] <- 0
  p <- egcar_prepare(X)
  for (penalty in c("l11", "l21")) for (lambda in c(0, 0.03, 10)) {
    fits <- lapply(c("cpp", "R", "reference"), function(backend) egcar_fit(
      p, rank = 1, penalty = penalty, lambda = lambda,
      control = egcar_control(backend = backend, max_iter = 35,
        abs_tol = 0, rel_tol = 0, adaptive_mu = TRUE, check_every = 5)))
    close(unlist(fits[[1]]$C), unlist(fits[[2]]$C))
    close(unlist(fits[[1]]$C), unlist(fits[[3]]$C))
    close(fits[[1]]$solver$primal_residual, fits[[3]]$solver$primal_residual)
    close(fits[[1]]$solver$dual_residual, fits[[3]]$solver$dual_residual)
    stopifnot(fits[[1]]$iterations == fits[[3]]$iterations)
    stopifnot(if (penalty == "l11") fits[[1]]$lambda_g == 0 else fits[[1]]$rho_e == 0)
    if (penalty == "l21") stopifnot(is.null(fits[[1]]$solver$Z), is.null(fits[[1]]$solver$H))
    warm <- lapply(c("cpp", "R", "reference"), function(backend) egcar_fit(
      p, rank = 1, penalty = penalty, lambda = 0.02,
      init = fits[[1]], control = egcar_control(backend = backend, max_iter = 10,
        abs_tol = 0, rel_tol = 0, adaptive_mu = FALSE)))
    close(unlist(warm[[1]]$C), unlist(warm[[2]]$C))
    close(unlist(warm[[1]]$C), unlist(warm[[3]]$C))
  }
}

fit <- egcar_fit(prepared, 2, "l21", lambda = 0,
                 control = egcar_control(max_iter = 1000))
stopifnot(!is.null(fit$L), identical(dim(coef(fit)), c(15L, 2L)))
stopifnot(length(coef(fit, by_view = TRUE)) == 3L)
stopifnot(length(predict(fit, views)) == 3L)
close(predict(fit, views, "sum"), Reduce(`+`, predict(fit, views)))
close(predict(fit, views[c(3, 1, 2)], "sum"), predict(fit, views, "sum"))
assert_error(egcar_fit(prepared, 1, "l21", init = list(C = list(matrix(0, 1, 1)))))
assert_error(egcar_fit(prepared, 1, "l11", init = fit))

# Rates retain the supplied script's formulas, and each uses only one penalty.
for (penalty in c("l11", "l21")) {
  rate <- egcar_rate(prepared, 1, penalty, multiplier = 2,
                     control = egcar_control(max_iter = 20))
  expected <- 2 * if (penalty == "l11") sqrt(log(15)/35) else sqrt((11 + log(15))/35)
  close(rate$lambda, expected)
}

# Shared folds, training-mean centering, full-grid scoring and loss direction.
shared <- egcar_cv_data(views, nfolds = 3, seed = 16)
for (fo in shared$folds) {
  close(fo$train_means[[1]], colMeans(views[[1]][fo$train_idx, , drop = FALSE]))
  close(colMeans(fo$train_views[[1]]), rep(0, 4), 1e-12)
}
before <- .Random.seed
cv <- egcar_cv(shared, 1, "l11", lambda = c(0.05, 0.01),
                control = egcar_control(max_iter_cv = 300, max_iter = 600))
stopifnot(identical(before, .Random.seed), identical(cv$fold_id, shared$fold_id))
stopifnot(nrow(cv$cv_fold_table) == 6, nrow(cv$cv_table) == 2)
close(cv$cv_table$mean_loss, -cv$cv_table$mean_score)
stopifnot(!is.null(cv$best), sum(cv$cv_table$selected) == 1)
close(cv$best$mean_loss, min(cv$cv_table$mean_loss))
stopifnot(is.finite(egcar_score(cv, views)))

# All-candidate failure is explicit, not replaced by an arbitrary penalty.
zero <- list(matrix(0, 8, 3), matrix(0, 8, 4))
failed <- egcar_cv(zero, 1, "l21", lambda = c(0.1, 0.01), nfolds = 2)
stopifnot(identical(failed$status, "no_valid_cv"), is.null(failed$fit), is.null(failed$best))
assert_error(coef(failed))
stopifnot(all(failed$cv_table$valid_folds == 0))

# Oracles remain separate truth-informed diagnostics.
oracle2 <- egcar_oracle_support(prepared, sim$active_local, 2)
stopifnot(!is.null(oracle2$L), identical(oracle2$penalty, "oracle-support"))
oracle1 <- egcar_oracle_population(sim$population, 2)
stopifnot(oracle1$rho_e == 0, oracle1$lambda_g == 0, !is.null(oracle1$solver$Z))
stopifnot(is.finite(egcar_distance(oracle1$L, sim$truth, 2)))
cat("Core package checks passed.\n")

# Positive-only EGCAR CV, but zero-coefficient fixed fits remain defined.
assert_error(EGCAR_L11_CV(shared, lambda = c(0.01, 0)))
assert_error(EGCAR_L21_CV(shared, lambda = c(0.01, 0)))
assert_error(egcar_experiment_config(rho_e_cv_grid = c(0.01, 0)))
assert_error(egcar_experiment_config(lambda_g_cv_grid = c(0.01, 0)))
stopifnot(all(egcar_experiment_config()$rho_e_cv_grid > 0),
          all(egcar_experiment_config()$lambda_g_cv_grid > 0))
