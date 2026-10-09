# Ten-method illustration on ONE dataset, using the same CV folds.
# Install the optional solvers first: install.packages(c("RGCCA", "PMA"))
# Then: source(system.file("examples", "compare_methods.R", package = "egcar"))
# This small example is not the full cluster benchmark.
library(egcar)
required <- c("RGCCA", "PMA")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install optional packages first: ", paste(missing, collapse = ", "))

sim <- egcar_simulate(n = 60, p_list = c(4, 4, 4), rank = 1,
                      active_per_view = 2, signal = 0.8, seed = 301)
shared <- egcar_cv_data(sim$views, nfolds = 3, seed = 302)
lambda_grid <- c(0.001, 0.01, 0.1)
workers <- 1L  # Fold-level parallelism: set to 3 after installing future/future.apply.

# Preserve SGCA's initializer/TGD defaults; shorten only the demonstration's
# WHOLE-CV time budget. A timeout remains a timeout, never a successful fit.
benchmarks <- benchmark_control(sgca_time_limit = 120)

fits <- list(
  Oracle1_population = egcar_oracle_population(sim$population, rank = 1),
  Oracle2_support = egcar_oracle_support(shared, sim$active_local, rank = 1),
  EGCAR_L11_rate = EGCAR_L11_Rate(shared, rank = 1),
  EGCAR_L11_CV = EGCAR_L11_CV(shared, rank = 1, lambda = lambda_grid,
                             workers = workers),
  EGCAR_L21_rate = EGCAR_L21_Rate(shared, rank = 1),
  EGCAR_L21_CV = EGCAR_L21_CV(shared, rank = 1, lambda = lambda_grid,
                             workers = workers),
  RGCCA = rgcca_cv(shared, rank = 1, tau_grid = c(0.1, 1),
                   workers = workers, benchmarks = benchmarks),
  SGCCA = sgcca_cv(shared, rank = 1, sparsity_grid = c(0.75, 1),
                   workers = workers, benchmarks = benchmarks),
  MultiCCA = multicca_cv(shared, rank = 1, penalty_grid = c(1.5, 2),
                         workers = workers, benchmarks = benchmarks),
  SGCA = sgca_cv(shared, rank = 1, k_grid = c(5, 10),
                 workers = workers, benchmarks = benchmarks)
)

report <- data.frame(
  method = names(fits),
  status = vapply(fits, function(f) f$status, character(1)),
  converged = vapply(fits, function(f) isTRUE(f$converged), logical(1)),
  has_loadings = vapply(fits, function(f) !is.null(f$L), logical(1)),
  fit_seconds = vapply(fits, function(f) f$fit_time, numeric(1)),
  total_seconds = vapply(fits, function(f) f$total_time, numeric(1)),
  row.names = NULL
)
print(report, row.names = FALSE)

# Oracles use population truth or true supports, unavailable in real data.
# SGCA uses its paper score by default; the other CV wrappers use the common
# generalized Rayleigh score. Same folds do not imply matched accuracy.
# For SGCA's common-score variant, explicitly set
# benchmark_control(sgca_cv_score = "common_loss", sgca_time_limit = 120).
# Methods run sequentially here; workers parallelizes folds within each CV.
