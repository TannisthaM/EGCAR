# Small EGCAR example. After installing the package, run:
# source(system.file("examples", "quickstart.R", package = "egcar"))
library(egcar)

# Matched observations in three views; one shared latent direction.
sim <- egcar_simulate(
  n = 90, p_list = c(6, 6, 6), rank = 1, active_per_view = 3,
  toeplitz_rho = c(0.5, 0.7, 0.9), signal = 0.8, seed = 101
)
views <- sim$views

# Prepare once so all methods use the same data and CV folds.
shared <- egcar_cv_data(views, nfolds = 3, seed = 102)
lambda_grid <- c(0.001, 0.01, 0.1)  # Small demonstration grid.
workers <- 1L  # Set to 3 with future and future.apply installed.

fits <- list(
  L11_fixed = EGCAR_L11(shared, rank = 1, lambda = 0.01),
  L21_fixed = EGCAR_L21(shared, rank = 1, lambda = 0.01),
  L11_rate = EGCAR_L11_Rate(shared, rank = 1),
  L21_rate = EGCAR_L21_Rate(shared, rank = 1),
  L11_CV = EGCAR_L11_CV(shared, rank = 1, lambda = lambda_grid,
                       workers = workers),
  L21_CV = EGCAR_L21_CV(shared, rank = 1, lambda = lambda_grid,
                       workers = workers)
)

# Keep completion and numerical convergence visible for every method.
summary_table <- data.frame(
  method = names(fits),
  status = vapply(fits, function(f) f$status, character(1)),
  converged = vapply(fits, function(f) isTRUE(f$converged), logical(1)),
  fit_seconds = vapply(fits, function(f) f$fit_time, numeric(1)),
  total_seconds = vapply(fits, function(f) f$total_time, numeric(1)),
  row.names = NULL
)
print(summary_table, row.names = FALSE)

# Inspect the CV-selected L21 result before extracting its loadings.
fit <- fits$L21_CV
print(fit)
if (identical(fit$status, "ok") && isTRUE(fit$converged)) {
  loadings <- coef(fit, by_view = TRUE)  # One p_k x rank matrix per view.
  scores <- predict(fit, newdata = views, type = "views")
  print(lapply(loadings, dim))
  print(head(scores[[1]]))  # In-sample scores, not a test-set evaluation.
}
# fit_time already includes loading extraction. For CV, total_time also
# includes tuning. Preparation of 'shared' is reported separately.
