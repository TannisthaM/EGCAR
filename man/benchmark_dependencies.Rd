\name{benchmark_dependencies}
\alias{benchmark_dependencies}
\title{Report installed optional comparison dependencies}
\description{Check installed package versions and the expected nonexported PMA helper names without installing packages.}
\usage{
benchmark_dependencies()
}
\value{A data frame of packages, installation status, versions, expected-helper presence and purpose. Helper presence is not a guarantee of numerical compatibility.}
\examples{
benchmark_dependencies()
}
