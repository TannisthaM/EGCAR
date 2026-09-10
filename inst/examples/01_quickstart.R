library(egcar)
sim <- egcar_simulate(n = 60, p_list = c(6, 7, 8), rank = 2,
                      active_per_view = 3, seed = 101)
shared <- egcar_cv_data(sim$views, nfolds = 3, seed = 102)
ctl <- egcar_control(backend = "cpp")
positive_grid <- c(0.001, 0.01, 0.1)
fits <- list(
  EGCAR_L11_CV = EGCAR_L11_CV(shared, rank = 2, lambda = positive_grid, control = ctl),
  EGCAR_L21_CV = EGCAR_L21_CV(shared, rank = 2, lambda = positive_grid, control = ctl),
  EGCAR_L11_Rate = EGCAR_L11_Rate(shared, rank = 2, control = ctl),
  EGCAR_L21_Rate = EGCAR_L21_Rate(shared, rank = 2, control = ctl))
lapply(fits, print)
# coef(fits$EGCAR_L21_CV, by_view = TRUE) returns one loading matrix per view,
# provided the CV fit has a valid requested-rank loading matrix.
