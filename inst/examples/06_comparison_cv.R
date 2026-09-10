library(egcar)
sim <- egcar_simulate(n = 60, p_list = c(8, 8, 8), rank = 1,
                      active_per_view = 3, seed = 301)
shared <- egcar_cv_data(sim$views, nfolds = 3, seed = 302)
bc <- benchmark_control()
# Same fold objects and training-mean centering in every call.
sgca <- sgca_cv(shared, rank = 1, k_grid = c(6, 12, 24),
  rho_grid = c(0, 0.01, 0.1), lambda_grid = c(0.01, 0.1, 1), benchmarks = bc)
rgcca <- rgcca_cv(shared, rank = 1, tau_grid = c(0.1, 0.5, 1), benchmarks = bc)
sgcca <- sgcca_cv(shared, rank = 1, sparsity_grid = c(0.5, 0.75, 1), benchmarks = bc)
multicca <- multicca_cv(shared, rank = 1, penalty_grid = c(1.5, 2, sqrt(8)), benchmarks = bc)
lapply(list(SGCA = sgca, RGCCA = rgcca, SGCCA = sgcca, MultiCCA = multicca), print)
