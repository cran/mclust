\name{one.XXX}
\alias{one.XI}

\title{Log-likelihood for a single cluster}
\description{
  \code{one.XXX} calculates the log-likelihood for a single
  unconstrained Gaussian cluster; \code{one.XI} calculates the
  log-likelihood for a single Gaussian spherical cluster.
}
\usage{
one.XXX(data)
one.XI(data)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{}
}

\examples{
data(iris)
one.XI(iris[,1:4])
one.XXX(iris[,1:4])
}
\keyword{clustering}
