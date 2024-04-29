\name{softmax}
\alias{softmax}

\title{Softmax function}

\description{
Efficient implementation (via Fortran) of the softmax (aka multinomial logistic) function converting a set of numerical values to probabilities summing to 1.
}

\usage{
softmax(x, v = NULL)
}

\arguments{
\item{x}{a matrix of dimension \eqn{n \times k} of numerical values. If a vector is provided, it is converted to a single-row matrix.}
\item{v}{an optional vector of length \eqn{k} of numerical values to be added to each row of \code{x} matrix. If not provided, a vector of zeros is used.}
}

\details{
Given the matrix \code{x}, for each row \eqn{x_{[i]} = [x_1, \dots, x_k]} (with \eqn{i=1,\dots,n}), the softmax function calculates
\deqn{
\text{softmax}(x_{[i]})_j = 
\dfrac{\exp(x_j + v_j)}{\sum_{l=1}^k \exp(x_l + v_l)}
\qquad \text{for } j = 1,\dots,k 
}
}

\value{Returns a matrix of the same dimension as \code{x} with values in the range \eqn{(0,1)} that sum to 1 along the rows.}

\author{Luca Scrucca}

\seealso{\code{\link{logsumexp}}}

\references{
Blanchard P., Higham D. J., Higham N. J. (2021). 
Accurately computing the log-sum-exp and softmax functions.
\emph{IMA Journal of Numerical Analysis}, 41/4:2311–2330.
\doi{10.1093/imanum/draa038}
}

\examples{
x = matrix(rnorm(15), 5, 3)
v = log(c(0.5, 0.3, 0.2))
(z = softmax(x, v))
rowSums(z)
}
