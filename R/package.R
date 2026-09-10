# Package-level documentation and explicit namespace declarations.
# No simulations, file writes, package installation or future plans run on load.
#' Accelerated sparse generalized correlation analysis
#'
#' @description Compiled and R implementations of separate entrywise and
#'   view-wise row-group regression estimators, with reusable common-loss
#'   cross-validation and optional comparison wrappers.
#' @details Start with \code{egcar_cv_data}, \code{egcar_fit},
#'   \code{egcar_rate} and \code{egcar_cv}. Select \code{penalty = "l11"}
#'   or \code{penalty = "l21"}; these penalties are not combined.
#'   The compiled backend is built at installation, not at fit time.
#'   Source organization follows the method-specific structure of ccar3;
#'   ccar3 is not a runtime dependency.
#' @name egcar-package
#' @aliases egcar
#' @docType package
#' @keywords internal
#' @useDynLib egcar, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict rnorm setNames
#' @importFrom utils getFromNamespace head tail
NULL
