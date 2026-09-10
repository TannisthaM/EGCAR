\name{egcar-package}
\alias{egcar}
\alias{egcar-package}
\title{Accelerated sparse generalized correlation analysis}
\description{Compiled and R implementations of separate entrywise and view-wise row-group regression estimators, with reusable common-loss CV and optional comparison wrappers.}
\details{Start with \code{egcar_simulate}, \code{egcar_cv_data},
\code{egcar_fit}, \code{egcar_rate} and \code{egcar_cv}.
The package compiles its C++14 backend during installation; it never calls
\code{sourceCpp} or runs a simulation on load. It does not depend on ccar3 at
runtime. Computational provenance, limitations, installation instructions and
all-method examples are installed under \code{doc} and \code{examples}.
The public names distinguish l11-only and l21-only estimators; old
The public interface exposes separate L11-only and L21-only families.}
