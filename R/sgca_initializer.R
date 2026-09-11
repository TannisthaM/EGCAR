# Required initializer dependency closure, bundled from:
# https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R
# Inspected 2026-09-10. Function bodies below are retained from that source.
# These four internal functions do not generate folds or perform CV or TGD.
# The experiment's common-loss CV and penalized TGD remain authoritative.
# Copyright and permission notice for this bundled code:
# MIT License
# 
# Copyright (c) 2026 Claire Donnat
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  sign(a)*pmax(0,abs(a)-b)
}

updatePi <- function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
  C <- Pi + 1/tau*A - nu/tau*B%*%Pi%*%B + nu/tau*sqB%*%(H-Gamma)%*%sqB
  D <- rho/tau
  Soft(C,D)
}
updateH <- function(sqB,Gamma,nu,Pi,K){
  temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
  temp <- (temp+t(temp))/2
  ev <- eigen(temp, symmetric = TRUE)
  d <- ev$values

  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(ev$vectors%*%diag(dfinal)%*%t(ev$vectors))
  }
  fr <- function(x) sum(pmin(1,pmax(d-x,0)))
  knots <- unique(c((d-1), d))
  knots <- sort(knots, decreasing=TRUE)
  temp2 <- which(sapply(knots, fr) <= K)
  lentemp <- tail(temp2, 1)
  a <- knots[lentemp]
  b <- knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a + (b-a) * (K-fa)/(fb-fa)
  dfinal <- pmin(1,pmax(d-theta,0))
  ev$vectors%*%diag(dfinal)%*%t(ev$vectors)
}
sgca_init_fixed <- function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
  A <- (A + t(A))/2
  B <- (B + t(B))/2
  p <- nrow(B)

  evB <- eigen(B, symmetric = TRUE)
  vals <- pmax(evB$values, 0)
  sqB  <- evB$vectors %*% diag(sqrt(vals), p, p) %*% t(evB$vectors)

  tau <- 4 * nu * (max(vals)^2)
  if (!is.finite(tau) || tau <= 0) tau <- 1

  criteria <- Inf
  iter <- 0
  H <- Pi <- oldPi <- diag(1, p)
  Gamma <- matrix(0, p, p)

  while(criteria > epsilon && iter < maxiter){
    for (j in 1:20){
      Pi <- updatePi(B, sqB, A, H, Gamma, nu, rho, Pi, tau)
    }
    H <- updateH(sqB, Gamma, nu, Pi, K)
    Gamma <- Gamma + (sqB %*% Pi %*% sqB - H) * nu
    criteria <- sqrt(sum((Pi - oldPi)^2))
    oldPi <- Pi
    iter <- iter + 1
    if (trace) cat("iter:", iter, "crit:", criteria, "\n")
  }
  list(Pi=Pi,H=H,Gamma=Gamma,iteration=iter,convergence=criteria)
}
