test_that("global row solver matches explicit ADMM and a closed-form solution", {
ns<-asNamespace('egcar')
eng<-function(backend) get('.egcar_engine',ns)(egcar_control(backend=backend))
a<-eng('cpp'); b<-eng('R')

set.seed(99)
for(dims in list(c(4,5,3),c(15,12,14),rep(8,5))) {
 n<-if(max(dims)>10) 7L else 40L
 views<-lapply(dims,function(p) matrix(rnorm(n*p),n,p))
 prep<-a$prepare_problem(a$center_views(views)$views)
 for(adapt in c(FALSE,TRUE)) {
  args<-list(prep=prep,lambda_g=.05,max_iter=250L,abs_tol=1e-8,rel_tol=1e-7,
    adaptive_mu=adapt,check_every=5L,keep_history=TRUE,compact_state=TRUE)
  x<-do.call(a$fit_l21_admm,args); y<-do.call(b$fit_l21_admm,args)
  stopifnot(max(abs(unlist(x$C_hat)-unlist(y$C_hat)))<1e-9,
    max(abs(as.matrix(x$history)-as.matrix(y$history)))<1e-7,
    x$iterations==y$iterations)
  full<-a$assemble_full_C(x$C_hat,prep$p_list)
  stopifnot(identical(full,t(full)))
  init<-a$egcar_warm_state(x,TRUE)
  args$init<-init; args$lambda_g<-.03;args$max_iter<-20L
  xx<-do.call(a$fit_l21_admm,args);yy<-do.call(b$fit_l21_admm,args)
  stopifnot(max(abs(unlist(xx$C_hat)-unlist(yy$C_hat)))<1e-9)
 }
 cat('native/R + history + adaptive mu + warm start:',paste(dims,collapse=','),'n=',n,'PASS\n')
}
# Identity block Hessians yield an independent closed-form answer:
# every directed row is soft group-thresholded at 2*lambda.
prep<-a$prepare_problem(lapply(1:3,function(k) matrix(rnorm(50*3),50,3)))
prep$S_kk<-rep(list(diag(3)),3);prep$eig<-rep(list(list(values=rep(1,3),vectors=diag(3))),3)
prep$S_kl<-lapply(prep$S_kl,function(S) S*3)
prep$egcar_context<-NULL;prep<-a$cache_problem_matrices(prep)
S<-a$assemble_full_C(prep$S_kl,prep$p_list);lam<-.075
Z<-S*pmax(0,1-2*lam/pmax(sqrt(rowSums(S*S)),.Machine$double.eps))
expected<-(Z+t(Z))/2
x<-a$fit_l21_admm(prep,lam,max_iter=3000L,abs_tol=1e-11,rel_tol=1e-10,
                  group_zero_tol=0,entry_zero_tol=0,keep_state=TRUE,compact_state=TRUE)
actual<-a$assemble_full_C(x$C_hat,prep$p_list)
stopifnot(x$converged,max(abs(actual-expected))<1e-8)
cat('closed-form global row penalty + final symmetry: PASS\n')
# L11 optimized path is unaffected by the new optional native branch.
old<-get('fit_l11_admm_reference',ns)(prep,.01,max_iter=100L)
new<-a$fit_l11_admm(prep,.01,max_iter=100L)
stopifnot(max(abs(unlist(old$C_hat)-unlist(new$C_hat)))<1e-9)
cat('L11 independent reference equivalence: PASS\n')
cat('ALL DIRECT SOLVER CHECKS PASSED\n')

  expect_true(x$converged)
})
