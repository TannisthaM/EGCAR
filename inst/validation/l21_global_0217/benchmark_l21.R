# Run against an installed egcar 0.2.17; arguments: old|new p_per_view blocks n.
library(egcar);ns<-asNamespace('egcar')
a<-commandArgs(TRUE);mode<-a[1];p<-as.integer(a[2]);K<-as.integer(a[3]);n<-as.integer(a[4])
e<-get('.egcar_engine',ns)(egcar_control())
old<-get('fit_l21_admm_legacy',ns);environment(old)<-e
set.seed(732);views<-lapply(seq_len(K),function(k) matrix(rnorm(n*p),n,p))
prep<-e$prepare_problem(e$center_views(views)$views)
f<-if(mode=='old') old else e$fit_l21_admm
lambda<-.05*sqrt(((K-1)*p+log(K*p))/n)
elapsed<-numeric(3);it<-numeric(3)
for(j in 1:3) {
 gc();t<-proc.time()[[3]]
 fit<-f(prep,lambda,max_iter=2000L,check_every=5L,keep_state=FALSE,compact_state=TRUE)
 load<-try(e$egcar_loading_from_operator(prep,fit$C_hat,1L,keep_full_C=FALSE),silent=TRUE)
 elapsed[j]<-proc.time()[[3]]-t;it[j]<-fit$iterations
 stopifnot(!inherits(load,"try-error"),isTRUE(load$valid))
}
cat(mode,p,K,n,'median_seconds',median(elapsed),'iterations',median(it),'converged',fit$converged,'\n')
