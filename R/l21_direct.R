# The fast native path stores right-rotated blocks and an exact scalar
# complement per global row. No coefficient or covariance direction is screened.
fit_l21_admm <- function(prep, lambda_g, mu_g = 1, max_iter = 2000L,
    abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
    scale_factor = 2, adapt_every = 10L, group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10, init = NULL, keep_history = FALSE,
    verbose = FALSE, check_every = 1L, keep_state = TRUE, compact_state = FALSE) {
  args <- as.list(environment())
  if(inherits(init,"egcar_fit")) init <- init$solver
  init <- .egcar_validate_init(init,prep,"l21")
  z <- egcar_get_context(prep)
  rotated <- is.null(init) || identical(init$state_basis,"right_spectral")
  if(EGCAR_BACKEND != "cpp" || !rotated) {
    args$init <- if(rotated && !is.null(init)) egcar_l21_expand_rotated(init,z) else init
    return(do.call(fit_l21_full_admm,args))
  }
  ctl <- egcar_controls(lambda_g,mu_g,max_iter,abs_tol,rel_tol,adaptive_mu,
    balance_ratio,scale_factor,adapt_every,check_every,keep_history)
  ctl$direct_l21 <- TRUE;ctl$rotated_l21 <- TRUE
  ctl$group_zero_tol <- group_zero_tol;ctl$entry_zero_tol <- entry_zero_tol
  ctl$keep_state <- keep_state
  ctl$zero_certificate <- lambda_g >= .5*max(unlist(egcar_group_norms(z,z$S),use.names=FALSE))*(1+16*.Machine$double.eps)
  E <- seq_along(z$edge_k)
  state <- if(is.null(init)) list(C=list(),cold=TRUE,Hk=list(),Hl=list(),
    a=lapply(z$p_list,function(p) rep(1,p)),hnull=lapply(z$p_list,numeric)) else
    c(list(C=list()),init[c("Hk","Hl","a","hnull")])
  raw <- egcar_native_solve(z,state,ctl,TRUE,verbose)
  C_hat <- raw$state$C;names(C_hat)<-z$keys
  active <- lapply(egcar_group_norms(z,C_hat),function(x) x>group_zero_tol)
  history <- if(keep_history) {
    h<-as.data.frame(raw$history);names(h)<-c("iter","objective","primal","dual","eps_primal","eps_dual","mu_g")
    h$iter<-as.integer(h$iter);h
  } else NULL
  out <- list(C_hat=C_hat,active_rows=active,
    directed_active_rows=lapply(raw$directed_norms,function(x) as.numeric(x)>group_zero_tol),
    formulation="global_rows_posthoc",state_basis="right_spectral",
    zero_certificate=ctl$zero_certificate,converged=raw$converged,iterations=raw$iterations,
    primal_residual=raw$primal,dual_residual=raw$dual,eps_primal=raw$eps_primal,
    eps_dual=raw$eps_dual,rho_e=0,lambda_g=lambda_g,mu_g=raw$mu,history=history,
    symmetrization_time=raw$symmetrization_time,
    convergence_scope="directed ADMM before final symmetrization")
  if(keep_state) {
    out$C<-C_hat
    out$Hk<-raw$state$Hk;out$Hl<-raw$state$Hl;out$a<-raw$state$a;out$hnull<-raw$state$hnull
  }
  out
}

# Diagnostic/reference conversion, never used inside the native iteration.
egcar_l21_expand_rotated <- function(init,z) {
  E<-seq_along(z$edge_k)
  hk<-lapply(E,function(e) {
    k<-z$edge_k[e];l<-z$edge_l[e];Q<-z$Q[[l]]
    A<-tcrossprod(init$Hk[[e]],Q)
    if(ncol(Q)<nrow(Q)) A<-A+(z$S[[e]]-tcrossprod(z$S[[e]]%*%Q,Q))*init$hnull[[k]]
    A
  })
  hl<-lapply(E,function(e) {
    k<-z$edge_k[e];l<-z$edge_l[e];Q<-z$Q[[k]]
    A<-Q%*%t(init$Hl[[e]])
    if(ncol(Q)<nrow(Q)) A<-A+sweep(z$S[[e]]-Q%*%crossprod(Q,z$S[[e]]),2,init$hnull[[l]],"*")
    A
  })
  list(Hk=hk,Hl=hl,a=init$a,formulation="global_rows_posthoc")
}

# Global row-L21 on a directed, zero-diagonal-block coefficient matrix.
# Minimize (1/4) sum_{k!=l} <C_kl, S_kk C_kl S_ll>
#          -(1/2) sum_{k!=l} <S_kl,C_kl> + lambda sum_i ||C_i,:||_2.
# The loss agrees with the former loss on symmetric C. Symmetry is imposed
# ONLY on the returned estimate by (Z+Z')/2 after ADMM. This projection is
# not asserted to minimize the symmetry-constrained objective.

fit_l21_full_admm <- function(prep, lambda_g, mu_g = 1, max_iter = 2000L,
    abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
    scale_factor = 2, adapt_every = 10L, group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10, init = NULL, keep_history = FALSE,
    verbose = FALSE, check_every = 1L, keep_state = TRUE, compact_state = FALSE) {
  ctl <- egcar_controls(lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu,
    balance_ratio, scale_factor, adapt_every, check_every, keep_history)
  ctl$direct_l21 <- TRUE
  z <- egcar_get_context(prep)
  if (!is.null(init) && !identical(init$formulation, "global_rows_posthoc"))
    stop("L21 0.2.17 warm starts must come from global_rows_posthoc; use init=NULL for old fits.")
  init <- .egcar_validate_init(init, prep, "l21")
  s <- if (is.null(init)) list(C=list(),
    Hk=empty_edge_list(prep$edge_table), Hl=empty_edge_list(prep$edge_table),
    a=lapply(z$p_list, function(pk) rep.int(1,pk))) else
    c(list(C=list()),init[c("Hk","Hl","a")])
  zero_bound <- .5*max(unlist(egcar_group_norms(z,z$S),use.names=FALSE))
  if(lambda_g >= zero_bound*(1+16*.Machine$double.eps)) {
    # Exact KKT certificate for zero: U=S/(2*mu), Z=0. This also supplies
    # a valid compact warm start without spending ADMM iterations on a
    # coefficient already known to have no loading directions.
    H <- lapply(z$S,function(A) A/(2*mu_g))
    ep <- sqrt(2*z$q)*abs_tol
    ed <- ep+rel_tol*sqrt(sum(vapply(z$S,function(A) sum(A*A),numeric(1)))/2)
    a <- list(state=list(Hk=H,Hl=H,a=lapply(z$p_list,numeric)),
      converged=TRUE,iterations=0L,primal=0,dual=0,eps_primal=ep,
      eps_dual=ed,mu=mu_g,history=if(keep_history)
        matrix(c(0,0,0,0,ep,ed,mu_g),1,7) else matrix(numeric(),0,7))
  } else a <- if (EGCAR_BACKEND == "cpp") egcar_native_solve(z,s,ctl,TRUE,verbose) else
    egcar_l21_direct_R(z,s,ctl,verbose)
  s <- a$state
  # The proximal variable Z carries row sparsity. Both directions are stored
  # in k-by-l orientation, so averaging never allocates the full P-by-P C.
  directed_norms <- Map(`*`,egcar_group_norms(z,s$Hk,s$Hl),s$a)
  directed_active <- lapply(directed_norms,function(x) x > group_zero_tol)
  t0 <- proc.time()[[3L]]
  C_hat <- lapply(seq_along(z$edge_k),function(e) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    ak <- s$a[[k]]*directed_active[[k]]; al <- s$a[[l]]*directed_active[[l]]
    A <- (s$Hk[[e]]*ak+sweep(s$Hl[[e]],2L,al,"*"))/2
    A[abs(A)<entry_zero_tol] <- 0
    A
  })
  names(C_hat) <- z$keys
  symmetry_time <- proc.time()[[3L]]-t0
  active <- lapply(egcar_group_norms(z,C_hat),function(x) x>group_zero_tol)
  h <- if (keep_history) {
    h <- as.data.frame(a$history)
    names(h) <- c("iter","objective","primal","dual","eps_primal","eps_dual","mu_g")
    h$iter <- as.integer(h$iter); h
  } else NULL
  out <- list(C_hat=C_hat,active_rows=active,directed_active_rows=directed_active,
    formulation="global_rows_posthoc",zero_certificate=(a$iterations==0L),converged=a$converged,iterations=a$iterations,
    primal_residual=a$primal,dual_residual=a$dual,eps_primal=a$eps_primal,
    eps_dual=a$eps_dual,rho_e=0,lambda_g=lambda_g,mu_g=a$mu,history=h,
    symmetrization_time=symmetry_time,
    convergence_scope="directed ADMM before final symmetrization")
  if (keep_state) {
    out$C <- C_hat # read-only alias of the final symmetric estimate, not another ADMM array
    out$Hk <- s$Hk;out$Hl <- s$Hl;out$a <- s$a
    # compact_state=FALSE no longer expands duplicate G/V; this model has a
    # single global proximal variable, represented by these row factors.
  }
  out
}

# Independent, explicit R implementation for checking the compressed native
# updates. It deliberately retains C/Z/U and is a diagnostic backend.
egcar_l21_direct_R <- function(z,s,ctl,verbose=FALSE) {
  E <- seq_along(z$edge_k); K <- seq_along(z$p_list)
  Zk <- Map(function(e) s$Hk[[e]]*s$a[[z$edge_k[[e]]]],E)
  Zl <- Map(function(e) sweep(s$Hl[[e]],2L,s$a[[z$edge_l[[e]]]],"*"),E)
  Uk <- Map(`-`,s$Hk,Zk); Ul <- Map(`-`,s$Hl,Zl)
  alpha <- 1.6
  mu <- ctl$mu; converged <- FALSE; history <- list()
  for(it in seq_len(ctl$max_iter)) {
    shift <- 2*mu
    update <- function(T,e) {
      Pt <- egcar_project(z,e,T)
      Ct <- (z$St[[e]]+shift*Pt)/(z$D[[e]]+shift)
      C <- if(z$full[[e]]) egcar_lift(z,e,Ct) else
        T+z$remainder[[e]]/shift+egcar_lift(z,e,Ct-Pt)
      list(C=C,Ct=Ct)
    }
    ck <- lapply(E,function(e) update(Zk[[e]]-Uk[[e]],e))
    cl <- lapply(E,function(e) update(Zl[[e]]-Ul[[e]],e))
    Ck <- lapply(ck,`[[`,"C"); Cl <- lapply(cl,`[[`,"C")
    Wk <- Map(function(C,Z,U) alpha*C+(1-alpha)*Z+U,Ck,Zk,Uk)
    Wl <- Map(function(C,Z,U) alpha*C+(1-alpha)*Z+U,Cl,Zl,Ul)
    norms <- egcar_group_norms(z,Wk,Wl)
    an <- lapply(norms,function(x) pmax(0,1-ctl$penalty/mu/pmax(x,.Machine$double.eps)))
    Nk <- lapply(E,function(e) Wk[[e]]*an[[z$edge_k[[e]]]])
    Nl <- lapply(E,function(e) sweep(Wl[[e]],2L,an[[z$edge_l[[e]]]],"*"))
    Uk <- Map(`-`,Wk,Nk); Ul <- Map(`-`,Wl,Nl)
    check <- it==1L || it==ctl$max_iter || it%%ctl$check_every==0L
    if(check) {
      ss <- function(x) sum(vapply(x,function(a) sum(a*a),numeric(1)))
      rp <- sqrt(ss(Map(`-`,Ck,Nk))+ss(Map(`-`,Cl,Nl)))
      rd <- mu*sqrt(ss(Map(`-`,Nk,Zk))+ss(Map(`-`,Nl,Zl)))
      ep <- sqrt(2*z$q)*ctl$abs_tol+ctl$rel_tol*max(sqrt(ss(Ck)+ss(Cl)),sqrt(ss(Nk)+ss(Nl)))
      ed <- sqrt(2*z$q)*ctl$abs_tol+ctl$rel_tol*mu*sqrt(ss(Uk)+ss(Ul))
      if(ctl$history) {
        obj <- sum(vapply(E,function(e) .25*sum(z$D[[e]]*(ck[[e]]$Ct^2+cl[[e]]$Ct^2))-
          .5*sum(z$S[[e]]*(Ck[[e]]+Cl[[e]])),numeric(1)))+
          ctl$penalty*sum(unlist(egcar_group_norms(z,Ck,Cl)))
        history[[length(history)+1L]] <- c(it,obj,rp,rd,ep,ed,mu)
      }
    }
    Zk <- Nk; Zl <- Nl
    if(check && rp<=ep && rd<=ed) {converged <- TRUE;break}
    if(check && ctl$adaptive && it%%ctl$adapt_every==0L) {
      nr<-rp/max(sqrt(ss(Ck)+ss(Cl)),sqrt(ss(Nk)+ss(Nl)),.Machine$double.eps)
      ns<-rd/max(mu*sqrt(ss(Uk)+ss(Ul)),.Machine$double.eps)
      factor <- if(nr>ctl$balance_ratio*max(ns,.Machine$double.eps)) ctl$scale_factor else
        if(ns>ctl$balance_ratio*max(nr,.Machine$double.eps)) 1/ctl$scale_factor else 1
      if(factor!=1) {mu<-mu*factor;Uk<-lapply(Uk,`/`,factor);Ul<-lapply(Ul,`/`,factor)}
    }
  }
  Hk <- Map(`+`,Zk,Uk);Hl <- Map(`+`,Zl,Ul)
  hn <- egcar_group_norms(z,Hk,Hl); zn <- egcar_group_norms(z,Zk,Zl)
  aa <- Map(function(a,b) ifelse(b>0,a/b,as.numeric(ctl$penalty==0)),zn,hn)
  list(state=list(Hk=Hk,Hl=Hl,a=aa),converged=converged,iterations=it,
    primal=rp,dual=rd,eps_primal=ep,eps_dual=ed,mu=mu,
    history=if(length(history)) do.call(rbind,history) else matrix(numeric(),0,7))
}
