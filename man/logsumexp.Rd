\name{logsumexp}
\alias{logsumexp}

\title{Log sum of exponentials}

\description{
Efficient implementation (via Fortran) of the log-sum-exp function.
}

\usage{
logsumexp(x, v = NULL)
}

\arguments{
\item{x}{a matrix of dimension \eqn{n \times k} of numerical values. If a vector is provided, it is converted to a single-row matrix.}
\item{v}{an optional vector of length \eqn{k} of numerical values to be added to each row of \code{x} matrix. If not provided, a vector of zeros is used.}
}

\details{
Given the matrix \code{x}, for each row \eqn{x_{[i]} = [x_1, \dots, x_k]} (with \eqn{i=1,\dots,n}), the log-sum-exp (LSE) function calculates
\deqn{
\text{LSE}(x_{[i]}) = \log \sum_{j=1}^k \exp(x_j + v_j) = m + \log \sum_{j=1}^k \exp(x_j + v_j - m)
}
where \eqn{m = \max(x_1+v_1, \dots, x_k+v_k)}.
}

\value{Returns a vector of values of length equal to the number of rows of \code{x}.}

\author{Luca Scrucca}

\seealso{\code{\link{softmax}}}

\references{
Blanchard P., Higham D. J., Higham N. J. (2021). 
Accurately computing the log-sum-exp and softmax functions.
\emph{IMA Journal of Numerical Analysis}, 41/4:2311–2330.
\doi{10.1093/imanum/draa038}
}

\examples{
x = matrix(rnorm(15), 5, 3)
v = log(c(0.5, 0.3, 0.2))
logsumexp(x, v)
}
