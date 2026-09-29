test_that("streamed edge reductions match full operators without input mutation", {
  set.seed(215)
  for (sizes in list(c(1L, 2L), c(3L, 1L, 5L), c(9L, 11L, 7L))) {
    tab <- egcar:::make_edge_table(sizes)
    z <- list(p_list = sizes, edge_k = tab$k, edge_l = tab$l)
    left <- lapply(seq_len(nrow(tab)), function(j) matrix(rnorm(tab$p_k[j]*tab$p_l[j]), tab$p_k[j]))
    names(left) <- tab$key
    right <- lapply(left, function(A) A * 2 + .1)
    before <- serialize(list(left, right), NULL)
    expected <- lapply(sizes, numeric)
    for(j in seq_along(left)) {
      expected[[tab$k[j]]] <- expected[[tab$k[j]]] + rowSums(left[[j]]^2)
      expected[[tab$l[j]]] <- expected[[tab$l[j]]] + colSums(right[[j]]^2)
    }
    expect_equal(egcar:::egcar_native_group_norms(left,right,tab$k,tab$l,sizes),
                 lapply(expected,sqrt),tolerance=1e-12)
    actual <- egcar:::egcar_native_edge_errors(left,right)
    A <- egcar:::assemble_full_C(left,sizes)
    B <- egcar:::assemble_full_C(right,sizes)
    expect_equal(unname(actual[c("error","reference")]),c(sqrt(sum((A-B)^2)),sqrt(sum(B^2))), tolerance=1e-12)
    expect_identical(serialize(list(left,right),NULL),before)
    expect_equal(egcar:::egcar_group_norms(z,left), egcar:::.egcar_engine(egcar_control(backend="R"))$egcar_group_norms(z,left),tolerance=1e-12)
  }
})

test_that("reduced population truth matches a dense eigenproblem", {
  for (K in 2:4) for (r in 1:3) {
    sim <- egcar_simulate(20,p_list=rep(5L,K),rank=r,active_per_view=3,
                          toeplitz_rho=.6,seed=100*K+r)
    pop <- sim$population
    half <- egcar:::matrix_power_psd(pop$Sigma0,.5)
    inverse <- egcar:::matrix_power_psd(pop$Sigma0,-.5)
    ev <- eigen(egcar:::symmetrize(half %*% pop$Cstar_full %*% half),symmetric=TRUE)
    dense <- inverse %*% ev$vectors[,seq_len(r),drop=FALSE]
    q1 <- egcar:::orthonormal_basis(dense);q2 <- egcar:::orthonormal_basis(sim$truth)
    expect_equal(pop$Sigma0_half,half,tolerance=1e-11)
    expect_equal(tcrossprod(q1),tcrossprod(q2),tolerance=1e-10)
    expect_equal(pop$eigenvalues,ev$values[seq_len(r)],tolerance=1e-10)
  }
})

test_that("fit timers include a deliberately delayed loading calculation", {
  sim <- egcar_simulate(35,p_list=c(4L,5L,6L),active_per_view=2,seed=215)
  prepared <- egcar_prepare(sim$views)
  ctl <- egcar_control()
  e <- egcar:::.egcar_engine(ctl)
  load <- e$egcar_loading_from_operator
  e$egcar_loading_from_operator <- function(...) { Sys.sleep(.06); load(...) }
  ans <- egcar:::.egcar_fit_prepared(prepared,1L,"l11",.01,ctl,e=e)
  expect_gte(ans$loading_time,.05)
  expect_gte(ans$fit_time,ans$solver_time+ans$loading_time-1e-6)
  expect_equal(ans$total_time,ans$fit_time)
  fit <- egcar_fit(prepared,lambda=.01)
  expect_gte(fit$wall_time,fit$fit_time)
  cv <- egcar_cv(sim$views,lambda=.01,nfolds=2)
  expect_equal(cv$total_time,cv$tuning_time+cv$fit_time,tolerance=1e-12)
  expect_gte(cv$fit_time,cv$loading_time+cv$solver_time-1e-6)
  expect_gte(cv$wall_time,cv$total_time-1e-6)
})

test_that("convergence eligibility is explicit and optional", {
  grid <- data.frame(candidate=1:2,lambda=c(.01,.1))
  rows <- data.frame(candidate=c(1,1,2,2),fold=c(1,2,1,2),
    loss=c(-10,-10,-2,-2),converged=c(TRUE,FALSE,TRUE,TRUE),iterations=1:4)
  strict <- egcar:::summarize_loading_cv(grid,rows,2,TRUE)
  legacy <- egcar:::summarize_loading_cv(grid,rows,2,FALSE)
  expect_equal(strict$eligible,c(FALSE,TRUE))
  expect_equal(strict$mean_loss,c(Inf,-2))
  expect_equal(legacy$mean_loss,c(-10,-2))
  expect_true(egcar_control()$cv_require_convergence)
  expect_true(egcar_experiment_config()$cv_require_convergence)
})

test_that("direct scoring reuses prediction and matches covariance scoring", {
  sim <- egcar_simulate(50,p_list=c(4,5,6),active_per_view=2,seed=216)
  fit <- egcar_fit(sim$views,lambda=.001)
  centered <- egcar:::center_views_at(sim$views,fit$means)
  val <- egcar:::make_validation_covariance(centered)
  expect_equal(egcar_score(fit,sim$views),egcar:::validation_score(fit$L,val),tolerance=1e-12)
  wide <- lapply(sim$views,function(X) X[1:3,,drop=FALSE])
  expect_equal(egcar_score(fit,wide),egcar:::validation_score(fit$L,
    egcar:::make_validation_covariance(egcar:::center_views_at(wide,fit$means))),tolerance=1e-12)
})

test_that("regularization paths retain numerical results while sharing factors", {
  sim <- egcar_simulate(35,p_list=c(4,5,6),active_per_view=2,seed=215)
  x <- egcar_prepare(sim$views)
  path <- egcar_path(x,lambda=c(.001,.01),penalty="l21")
  a <- egcar_fit(x,lambda=.01,penalty="l21")
  b <- egcar_fit(x,lambda=.001,penalty="l21",init=a)
  expect_equal(path$fits[[1]]$C,b$C,tolerance=1e-12)
  expect_equal(path$fits[[2]]$C,a$C,tolerance=1e-12)
  expect_equal(path$fits[[1]]$L,b$L,tolerance=1e-10)
})

test_that("experiment metrics agree with dense operator diagnostics", {
  sim <- egcar_simulate(25,p_list=c(4,5,6),active_per_view=2,seed=215)
  fit <- egcar_fit(sim$views,lambda=.01)
  row <- egcar:::evaluate_method("test",fit$C,fit$loading,sim$population,1,25,1,converged=fit$converged)
  C <- egcar:::assemble_full_C(fit$C,sim$population$p_list)
  expect_equal(row$C_error,egcar:::frob(C-sim$population$Cstar_full),tolerance=1e-11)
  sm <- egcar:::support_metrics(C,sim$population$active_global)
  expect_equal(row$support_precision,unname(sm[["precision"]]))
  expect_equal(row$support_recall,unname(sm[["recall"]]))
  expect_equal(row$support_fdp,unname(sm[["fdp"]]))
})

test_that("study CV timers include final loading and successful summaries filter failures", {
  sim <- egcar_simulate(35,p_list=c(4,5,6),active_per_view=2,seed=215)
  shared <- egcar_cv_data(sim$views,nfolds=2)
  e <- egcar:::.egcar_engine(egcar_control())
  for(nm in c('cross_validate_penalties','fit_estimator')) {
    fn <- getFromNamespace(nm,'egcar'); environment(fn)<-e; assign(nm,fn,e)
  }
  list2env(setNames(unclass(egcar_experiment_config()),toupper(names(egcar_experiment_config()))), e)
  load <- e$egcar_loading_from_operator
  e$egcar_loading_from_operator <- function(...) {Sys.sleep(.04);load(...)}
  for(method in c('l11','l21')) {
    ans <- e$cross_validate_penalties(sim$views,shared$full$prep,shared$folds,1,.01,.01,method)
    expect_identical(ans$status,'ok')
    expect_gte(ans$fit_time,.035)
  }
  rows <- data.frame(rank=1,n=35,method='EGCAR-L11-rate',status=c('ok','time_limit','not_converged'),total_time=c(1,99,999))
  summary <- egcar:::summarize_metric; environment(summary) <- e
  e$METHOD_ORDER <- unique(rows$method)
  expect_equal(summary(rows,'total_time',1)$mean,1)
})
