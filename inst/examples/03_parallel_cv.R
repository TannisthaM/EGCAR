library(egcar)
sim <- egcar_simulate(n = 60, p_list = c(6, 7, 8), active_per_view = 3, seed = 24)
shared <- egcar_cv_data(sim$views, nfolds = 3, seed = 25)
# All methods reuse this plan and the same folds. No method-level nested parallelism.
old <- future::plan()
tryCatch({
  future::plan(future::multisession, workers = 3L)
  a <- EGCAR_L11_CV(shared, lambda = c(0.001, 0.01, 0.1), workers = 3L)
  b <- EGCAR_L21_CV(shared, lambda = c(0.001, 0.01, 0.1), workers = 3L)
  print(a); print(b)
}, finally = future::plan(old))
