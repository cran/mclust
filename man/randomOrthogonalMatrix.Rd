\name{randomOrthogonalMatrix}
\alias{randomOrthogonalMatrix}

\title{Random orthogonal matrix}

\description{
Generate a random orthogonal basis matrix of dimension \eqn{(nrow x ncol)} using 
the method in Heiberger (1978).
}

\usage{
randomOrthogonalMatrix(nrow, ncol, n = nrow, d = ncol, seed = NULL)
}

\arguments{
  \item{nrow}{the number of rows of the resulting orthogonal matrix.}
  \item{ncol}{the number of columns of the resulting orthogonal matrix.}
  \item{n}{deprecated. See \code{nrow} above.}
  \item{d}{deprecated. See \code{ncol} above.}
  \item{seed}{an optional integer argument to use in \code{set.seed()} for 
              reproducibility. By default the current seed will be used.
              Reproducibility can also be achieved by calling \code{set.seed()}
              before calling this function.}
}

\details{
The use of arguments \code{n} and \code{d} is deprecated and they will be removed in the future.
}

\value{
An orthogonal matrix of dimension \eqn{nrow x ncol} such that each column is orthogonal to the other and has unit lenght. Because of the latter, it is also called orthonormal.
}

\seealso{\code{\link{coordProj}}}

\references{
Heiberger R. (1978) Generation of random orthogonal matrices. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 27(2), 199-206.
}

\examples{
B <- randomOrthogonalMatrix(10,3)
zapsmall(crossprod(B))
}
