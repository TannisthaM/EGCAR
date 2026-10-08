# Usage: Rscript benchmark_tier.R /path/to/Rlib p_per_view
a<-commandArgs(TRUE);stopifnot(length(a)==2L);p<-as.integer(a[2])
.libPaths(c(normalizePath(a[1]),.libPaths()))
library(egcar);mode<-as.character(packageVersion('egcar'));e<-get('.egcar_engine',asNamespace('egcar'))(egcar_control())
# Exact original tier 1/2/3 seed rules, rank 1, replication 1, signal 0.8.
r<-1L;rep<-1L;n<-150L;task<-20260907L+100000L*rep+1000L*r
popseed<-task+100000L+1000L*r;dataseed<-popseed+77L
pop<-e$make_population(rep(p,3L),r,5L,c(.5,.7,.9),.8,popseed)
v<-e$simulate_views(pop,n,dataseed);prep<-e$prepare_problem(e$center_views(v)$views)
cat('data_sha256',digest::digest(v,algo='sha256',serializeVersion=2),'\n')
set.seed(task+1000000L+10000L*r+n+5000L)
lambda<-.05*sqrt((2*p+log(3*p))/n)
gc();t<-proc.time()[3];fit<-e$fit_l21_admm(prep,lambda,max_iter=2000L,check_every=5L,keep_state=FALSE)
solver<-proc.time()[3]-t
load<-e$egcar_loading_from_operator(prep,fit$C_hat,r,keep_full_C=FALSE)
total<-proc.time()[3]-t
cat(mode,'p',p,'solver_seconds',solver,'fit_through_loadings_seconds',total,'iterations',fit$iterations,
 'converged',fit$converged,'loading_valid',load$valid,'\n')
