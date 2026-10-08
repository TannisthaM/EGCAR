test_that("rotated L21 preserves genuine covariance nullspace contributions", {
 ns<-asNamespace("egcar");eng<-function(backend)get(".egcar_engine",ns)(egcar_control(backend=backend))
 a<-eng("cpp");b<-eng("R");set.seed(288)
 prep<-a$prepare_problem(lapply(c(3L,4L,5L),function(p)matrix(rnorm(20*p),20,p)))
 prep$eig<-lapply(prep$p_list,function(p)list(values=c(1,rep(0,p-1)),vectors=diag(p)))
 prep$S_kk<-lapply(prep$p_list,function(p)diag(c(1,rep(0,p-1))))
 prep$S_kl<-lapply(prep$S_kl,function(M)matrix(rnorm(length(M)),nrow(M),ncol(M)))
 prep$egcar_context<-NULL;prep<-a$cache_problem_matrices(prep)
 for(lam in c(0,.03,2))for(adapt in c(FALSE,TRUE)) {
  args<-list(prep=prep,lambda_g=lam,max_iter=25L,abs_tol=0,rel_tol=0,
   adaptive_mu=adapt,check_every=1L,keep_history=TRUE,group_zero_tol=0,entry_zero_tol=0)
  x<-do.call(a$fit_l21_admm,args);y<-do.call(b$fit_l21_admm,args)
  expect_equal(unname(x$C_hat),unname(y$C_hat),tolerance=1e-9)
  expect_equal(as.matrix(x$history),as.matrix(y$history),tolerance=1e-8)
  expect_identical(x$state_basis,"right_spectral")
  expect_equal(vapply(x$Hk,ncol,1L),rep(1L,3))
  args$init<-a$egcar_warm_state(x,TRUE);args$max_iter<-10L
  xx<-do.call(a$fit_l21_admm,args);yy<-do.call(b$fit_l21_admm,args)
  expect_equal(unname(xx$C_hat),unname(yy$C_hat),tolerance=1e-8)
 }
})

test_that("zero certificates and compact states preserve the requested loading rank", {
 ns<-asNamespace("egcar");a<-get(".egcar_engine",ns)(egcar_control())
 set.seed(188);prep<-a$prepare_problem(lapply(1:3,function(k)matrix(rnorm(12*30),12,30)))
 x<-a$fit_l21_admm(prep,1e5,keep_state=TRUE)
 expect_true(x$zero_certificate);expect_true(x$converged);expect_equal(x$iterations,0L)
 expect_true(all(unlist(x$C_hat)==0));expect_identical(x$state_basis,"right_spectral")
 y<-a$fit_l21_admm(prep,.03,init=a$egcar_warm_state(x,TRUE),max_iter=1000L)
 expect_identical(dim(a$assemble_full_C(y$C_hat,prep$p_list)),c(90L,90L))
 z<-a$fit_l21_admm(prep,.03,keep_state=FALSE)
 expect_false(any(c("Hk","Hl","hnull","a") %in% names(z)))
 expect_true(z$symmetrization_time>=0)
})
