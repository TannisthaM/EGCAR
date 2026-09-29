# An independent, deliberately direct transcription of Algorithm 1.
paper_direct <- function(A, B, init, k, lambda, eta, steps) {
  hard <- function(Z) {
    j <- order(-rowSums(Z^2), seq_len(nrow(Z)))
    if (k < nrow(Z)) Z[j[-seq_len(k)], ] <- 0
    Z
  }
  power <- function(M, exponent) {
    e <- eigen((M + t(M))/2, symmetric=TRUE)
    e$vectors %*% diag(e$values^exponent, nrow(M)) %*% t(e$vectors)
  }
  U <- hard(init)
  U <- U %*% power(t(U) %*% B %*% U, -.5)
  V <- U %*% power(diag(ncol(U)) + t(U) %*% A %*% U / lambda, .5)
  for (i in seq_len(steps)) {
    previous <- V
    V <- hard(V - 2 * eta * (-A %*% V + lambda * B %*% V %*%
      (t(V) %*% B %*% V - diag(ncol(V)))))
  }
  list(L=V %*% power(t(V) %*% B %*% V, -.5), change=sqrt(sum((V-previous)^2)))
}

test_that("SGCA follows Algorithm 1 including rank-two thresholding and factor two", {
  set.seed(216)
  X <- matrix(rnorm(120), 20, 6)
  A <- crossprod(X)/20
  B <- diag(diag(A))
  for (r in 1:2) {
    U <- matrix(rnorm(6*r), 6, r)
    for (steps in c(1, 17, 63)) {
      expected <- paper_direct(A, B, U, 4, .01, .001, steps)
      actual <- sgca_tgd_penalized(A, B, U, r, 4, .01, max_iter=steps)
      expect_equal(actual$L, expected$L, tolerance=2e-10)
      expect_equal(actual$absolute_change, expected$change, tolerance=2e-10)
      expect_equal(crossprod(actual$L, B %*% actual$L), diag(r), tolerance=1e-10)
      expect_lte(sum(rowSums(actual$L^2)>0), 4)
      expect_equal(actual$iterations, steps)
    }
  }
  U <- matrix(1, 5, 2)
  expect_equal(sgca_hard_rows(U, 2), rbind(U[1:2,], matrix(0,3,2)))
})

test_that("fixed paper iterations do not stop early or assert convergence", {
  A <- diag(c(3,2,1)); B <- diag(3); U <- matrix(c(1,0,0),3,1)
  fixed <- sgca_tgd_penalized(A,B,U,1,3,.01)
  expect_equal(fixed$iterations,15000)
  expect_true(fixed$completed)
  expect_equal(fixed$stop_reason,"fixed_iterations")
  adaptive <- sgca_tgd_penalized(A,B,U,1,3,.01,max_iter=Inf,stopping="absolute_change")
  expect_equal(adaptive$iterations,1)
  expect_true(adaptive$converged)
  expect_error(benchmark_control(sgca_tgd_max_iter=Inf),"finite")
  expect_equal(benchmark_control(sgca_tgd_max_iter=Inf,sgca_stopping="absolute_change")$sgca_tgd_max_iter,Inf)
  expect_error(egcar_experiment_config(sgca_max_iter_tgd=Inf),"finite")
})

test_that("Fantope projection enforces trace equality even below the trace bound", {
  for (p in c(4L,9L)) for (r in c(1L,2L)) {
    H <- updateH(diag(p),matrix(0,p,p),1,diag(-2,p),r)
    e <- eigen(H,symmetric=TRUE)$values
    expect_equal(sum(e),r,tolerance=1e-10)
    expect_gte(min(e),-1e-10)
    expect_lte(max(e),1+1e-10)
    expect_equal(H,diag(r/p,p),tolerance=1e-10)
  }
  A <- diag(4)+.02; B <- diag(seq(1,2,length.out=4))
  for (nu in c(1,2)) {
    z <- sgca_init_fixed(A,B,.05,2,nu=nu,epsilon=1e-100,maxiter=7)
    cached <- sgca_init_cached(sgca_prepare_initializer(A,B,sgca_init_fixed,nu),.05,2,epsilon=1e-100,maxiter=7)
    expect_equal(cached,z,tolerance=1e-10)
    expect_equal(z$iteration,7)
    expect_equal(sum(diag(z$H)),2,tolerance=1e-10)
  }
})

test_that("paper CV recomputes rho per training sample and retains completion diagnostics", {
  sim <- egcar_simulate(n=31,p_list=c(3,3,3),active_per_view=2,seed=216)
  shared <- egcar_cv_data(sim$views,fold_id=rep(1:3,length.out=31))
  bc <- benchmark_control(sgca_time_limit=Inf,sgca_init_max_iter=2,sgca_tgd_max_iter=7)
  fit <- sgca_cv(shared,k_grid=c(5,9),benchmarks=bc)
  expect_equal(fit$status,"ok")
  expect_true(fit$completed)
  expect_false(fit$converged)
  expect_true(all(fit$cv_table$eligible))
  expect_true(all(fit$cv_fold_table$completed))
  expect_equal(fit$cv_fold_table$sgca_rho_actual,.5*sqrt(log(9)/fit$cv_fold_table$n_train))
  expect_equal(fit$diagnostics$sgca_rho_actual,.5*sqrt(log(9)/31))
  expect_true(all(fit$cv_fold_table$sgca_tgd_iterations==7))
  expect_equal(fit$diagnostics$sgca_init_stop_reason,"iteration_limit")
  expect_true(all(is.na(fit$cv_table$sgca_rho)))
  expect_true(all(fit$cv_table$sgca_lambda==.01))
  expect_equal(fit$total_time,fit$fit_time+fit$tuning_time)
  bc$sgca_cv_require_convergence <- TRUE
  strict <- sgca_cv(shared,k_grid=c(5,9),benchmarks=bc)
  expect_equal(strict$status,"no_valid_cv")
  expect_null(strict$L)
  bc$sgca_cv_require_convergence <- FALSE
  direct <- sgca_cv(shared,k_grid=5,rho_grid=.03,benchmarks=bc)
  expect_true(all(direct$cv_fold_table$sgca_rho_actual==.03))
  expect_equal(direct$diagnostics$sgca_rho_actual,.03)
})

test_that("paper validation score uses test covariance and is distinct from common loss", {
  set.seed(216)
  views <- list(matrix(rnorm(32)+5,8,4),matrix(rnorm(24)-3,8,3))
  L <- matrix(rnorm(14),7,2)
  validation <- make_validation_covariance(views)
  X <- scale(do.call(cbind,views),scale=FALSE)
  expected <- -sum(diag(t(L)%*%(crossprod(X)/8)%*%L))
  expect_equal(sgca_paper_validation_loss(L,validation),expected,tolerance=1e-12)
  expect_gt(abs(expected-validation_loss(L,validation)),.01)
  wide <- lapply(views, function(X) X[1:5,,drop=FALSE])
  val_wide <- make_validation_covariance(wide)
  X_wide <- scale(do.call(cbind,wide),scale=FALSE)
  expect_equal(sgca_paper_validation_loss(L,val_wide),
    -sum(diag(t(L)%*%(crossprod(X_wide)/5)%*%L)),tolerance=1e-12)
  expect_error(sgca_paper_validation_loss(L,list(Sigma=diag(7))),"recreate")
  sim <- egcar_simulate(n=24,p_list=c(3,3,3),active_per_view=2,seed=216)
  shared <- egcar_cv_data(sim$views,nfolds=2)
  b <- benchmark_control(sgca_time_limit=Inf,sgca_init_max_iter=5,sgca_tgd_max_iter=15,sgca_cv_score="common_loss")
  fit <- sgca_cv(shared,k_grid=6,benchmarks=b)
  expect_equal(fit$diagnostics$sgca_cv_score,"common_loss")
  expect_equal(fit$status,"ok")
})

test_that("paper defaults reach the experiment engine and retained compact outputs", {
  cfg <- egcar_experiment_config()
  bc <- benchmark_control()
  expect_null(cfg$sgca_rho_grid)
  expect_equal(cfg$sgca_k_grid,seq(5,100,5))
  expect_equal(cfg$sgca_lambda_grid,.01)
  expect_equal(bc$sgca_ridge,0)
  expect_equal(bc$sgca_init_max_iter,1000)
  expect_equal(bc$sgca_tgd_max_iter,15000)
  expect_equal(bc$sgca_stopping,"fixed_iterations")
  expect_equal(bc$sgca_cv_score,"paper")
  cfg <- egcar_experiment_config(p_list=c(3L,3L,3L),n_grid=18L,rank_grid=1L,
    active_per_view=2L,n_folds=2L,rho_e_cv_grid=.001,lambda_g_cv_grid=.001,
    max_iter_cv=5L,max_iter_final=5L,oracle1_max_iter=5L,
    sgca_k_grid=5L,sgca_max_iter_init=2L,sgca_max_iter_tgd=7L,sgca_time_limit=Inf,
    stop_if_benchmark_packages_missing=FALSE,make_plots=FALSE,make_loading_plots=FALSE,
    save_cv_fold_results=TRUE)
  outdir <- tempfile("sgca-paper-study-")
  on.exit(unlink(outdir,recursive=TRUE),add=TRUE)
  out <- suppressWarnings(run_egcar_experiments(outdir,config=cfg,backend="R"))
  row <- out$results[out$results$method=="SGCA",]
  expect_equal(row$status,"ok")
  expect_true(row$benchmark_completed)
  expect_equal(row$sgca_tgd_iterations,7)
  expect_equal(row$sgca_rho_actual,.5*sqrt(log(9)/18))
  expect_equal(out$config$sgca_stopping,"fixed_iterations")
  expect_equal(out$config$sgca_cv_score,"paper")
  saved <- read.csv(file.path(outdir,"simulation_results.csv"))
  expect_equal(saved$sgca_tgd_iterations[saved$method=="SGCA"],7)
})
