\name{covw}
\alias{covw}

\title{Weighted means, covariance and scattering matrices conditioning on a weighted matrix.}

\description{
Compute efficiently (via Fortran code) the means, covariance and scattering matrices conditioning on a weighted or indicator matrix
}
\usage{
covw(X, Z, normalize = TRUE)
}

\arguments{
  \item{X}{A \eqn{(n x p)} data matrix, with \eqn{n} observations on \eqn{p} variables.}
  \item{Z}{A \eqn{(n x G)} matrix of weights, with \eqn{G} number of groups.}
  \item{normalize}{A logical indicating if rows of \code{Z} should be normalized to sum to one.}
}

\value{A list with the following components:
  \item{mean}{A \eqn{(p x G)} matrix of weighted means.}
  \item{S}{A \eqn{(p x p x G)} array of weighted covariance matrices.}
  \item{W}{A \eqn{(p x p x G)} array of weighted scattering matrices.}
  
}

%\seealso{}
\author{M. Fop and L. Scrucca}

\examples{
# Z as an indicator matrix
X = iris[,1:4]
Z = unmap(iris$Species)
str(covw(X, Z))
# Z as a matrix of weights
mod = Mclust(X, G = 3, modelNames = "VVV")
str(covw(X, mod$z))
}

\keyword{multivariate}
