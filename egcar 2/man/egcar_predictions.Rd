\name{egcar_predictions}
\alias{egcar_score}
\alias{coef.egcar_fit}
\alias{coef.egcar_cv}
\alias{predict.egcar_fit}
\alias{predict.egcar_cv}
\title{Extract loadings, calculate scores and predict components}
\description{Access stacked or view-wise coefficients, compute component scores with stored training means, or evaluate the common generalized Rayleigh score.}
\usage{
egcar_score(object, newdata, ridge = 1e-8)
}
\arguments{
\item{object}{An \code{egcar_fit} or \code{egcar_cv} object with a valid loading matrix.}
\item{newdata}{Complete matched views with training feature names and column order. Named views may be supplied in a different view order and are matched by name.}
\item{ridge}{Relative score-normalization ridge, default 1e-8.}
}
\details{Use \code{coef(object)} for the stacked loading matrix or \code{coef(object, by_view = TRUE)} for its view blocks. Use \code{predict(object, newdata, type = "views")} for a list of per-view component scores, or \code{type = "sum"} to sum them. Prediction does not impute an absent view. \code{egcar_score} reports the positive score; CV minimizes its negative.}
\value{A scalar validation score. The coefficient and prediction methods return matrices or lists of matrices.}
