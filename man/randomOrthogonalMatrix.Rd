\name{randomOrthogonalMatrix}
\alias{randomOrthogonalMatrix}

\title{Random orthogonal matrix}

\description{
Generate a random orthogonal basis matrix of dimension \eqn{(n x d)} using 
the method in Heiberger (1978).
}

\usage{
randomOrthogonalMatrix(n, d)
}

\arguments{
  \item{n}{the number of rows of the resulting orthogonal matrix.}
  \item{d}{the number of columns of the resulting orthogonal matrix.}
}
\value{
An orthogonal matrix of dimension \eqn{n x d} such that each column is orthogonal to the other and has unit lenght.}

\seealso{\code{\link{coordProj}}}

\references{
Heiberger R. (1978) Generation of random orthogonal matrices. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 27(2), 199-206.
}

\examples{
B <- randomOrthogonalMatrix(10,3)
zapsmall(crossprod(B))
}
