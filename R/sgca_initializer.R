# Required initializer dependency closure, bundled from:
# https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R
# Inspected 2026-09-10; adapted in 0.2.16 to enforce the paper's trace equality.
# These four internal functions do not generate folds or perform CV or TGD.
# The paper-based SGCA CV and Algorithm 1 are in alt_SGCA.R.
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
  C <- Pi + 1/tau*A - nu/tau*B%*%Pi%*%B + nu/tau*sqB%*%(H-Gamma/nu)%*%sqB
  D <- rho/tau
  Soft(C,D)
}
updateH <- function(sqB,Gamma,nu,Pi,K){
  temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
  temp <- (temp+t(temp))/2
  ev <- eigen(temp, symmetric = TRUE)
  d <- ev$values

  # Projection onto {0 <= H <= I, trace(H) = K}, equation (14).
  # The authors' R helper instead accepts trace <= K when clipping suffices.
  if (K < 0 || K > length(d)) stop("Invalid Fantope trace.")
  if (K == 0) return(matrix(0, length(d), length(d)))
  if (K == length(d)) return(diag(length(d)))
  lower <- min(d) - 1
  upper <- max(d)
  for (j in seq_len(80L)) {
    theta <- (lower + upper) / 2
    clipped <- pmin(1, pmax(d - theta, 0))
    if (sum(clipped) > K) lower <- theta else upper <- theta
    if (abs(sum(clipped) - K) <= 1e-12 * max(1, K)) break
  }
  tcrossprod(sweep(ev$vectors, 2L, clipped, "*"), ev$vectors)
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
