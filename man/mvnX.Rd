\name{mvnX}
\alias{mvnX}
\alias{mvnXII}
\alias{mvnXXI}
\alias{mvnXXX}
\title{
  Multivariate Normal Fit
}
\description{
  Computes the mean, covariance, and loglikelihood from fitting a single
  MVN or Gaussian.
}
\usage{
mvnX(data)
mvnXII(data)
mvnXXI(data)
mvnXXX(data)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
}
\value{
  A list of including the parameters of the Gaussian model best fitting the
  data, and the corresponding loglikelihood for the data under the model.
}
\details{
  \code{mvnXII} computes the best fitting
  Gaussian with the covariance restricted to be a multiple of the identity. 
  \code{mvnXXI} computes the best fitting
  Gaussian with the covariance restricted to be diagonal. 
  \code{mvnXXX} computes the best fitting
  Gaussian with ellipsoidal (unrestricted) covariance. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{mvn}},
  \code{\link{mstepE}}
}
\examples{
n <- 1000

set.seed(0)
x <- rnorm(n, mean = -1, sd = 2)
mvnX(x) 

mu <- c(-1, 0, 1)

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% (2*diag(3)), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvnXII(x) 

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% diag(1:3), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvnXXI(x)

Sigma <- matrix(c(9,-4,1,-4,9,4,1,4,9), 3, 3)
set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% chol(Sigma), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvnXXX(x) 
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
