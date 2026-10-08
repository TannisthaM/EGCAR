# Usage: Rscript bench.R /path/to/Rlib p_per_view K n
a<-commandArgs(TRUE);stopifnot(length(a)==4L);p<-as.integer(a[2]);K<-as.integer(a[3]);n<-as.integer(a[4])
.libPaths(c(normalizePath(a[1]),.libPaths()))
library(egcar);mode<-as.character(packageVersion('egcar'));ns<-asNamespace('egcar');invisible(loadNamespace('RSpectra'));invisible(gc())
e<-get('.egcar_engine',ns)(egcar_control())

set.seed(732);views<-lapply(seq_len(K),function(k) matrix(rnorm(n*p),n,p))
prep<-e$prepare_problem(e$center_views(views)$views)
f<-e$fit_l21_admm
lambda<-.05*sqrt(((K-1)*p+log(K*p))/n)
elapsed<-numeric(3);it<-numeric(3);sol<-numeric(3);lt<-numeric(3)
for(j in 1:3) {
 gc();t<-proc.time()[[3]]
 fit<-f(prep,lambda,max_iter=2000L,check_every=5L,keep_state=FALSE,compact_state=TRUE)
 sol[j]<-proc.time()[[3]]-t; tload<-proc.time()[[3]]
 load<-try(e$egcar_loading_from_operator(prep,fit$C_hat,1L,keep_full_C=FALSE),silent=TRUE)
 elapsed[j]<-proc.time()[[3]]-t;lt[j]<-proc.time()[[3]]-tload;it[j]<-fit$iterations
 stopifnot(!inherits(load,"try-error"),isTRUE(load$valid))
}
cat(mode,p,K,n,'median_solver',median(sol),'median_loading',median(lt),'median_seconds',median(elapsed),'iterations',median(it),'converged',fit$converged,'\n')
