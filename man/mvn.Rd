\name{mvn}
\alias{mvn}
\title{
  Multivariate Normal Fit
}
\description{
  Computes the mean, covariance, and loglikelihood from fitting a single
  MVN or Gaussian to given data.
}
\usage{
mvn( modelName, data)
}
\arguments{
  \item{modelName}{
    A character string representing a model name. This can be either
    \code{"Spherical"}, \code{"Diagonal"}, or \code{"Ellipsoidal"} or an
    MCLUST-style model name: \cr 
    "E", "V", "X" (one-dimensional) \cr
    "EII", "VII", "XII" (spherical) \cr
    "EEI", "VEI", "EVI", "VVI", "XXI" (diagonal) \cr
    "EEE", "EEV", "VEV", "VVV", "XXX" (ellipsoidal) 
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
}
\value{
  A list of including the parameters of the Gaussian model best fitting the
  data, and the corresponding loglikelihood for the data under the model.
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
  \code{\link{mvnX}},
  \code{\link{mvnXII}},
  \code{\link{mvnXXI}},
  \code{\link{mvnXXX}},
  \code{\link{mstep}}
}
\examples{
n <- 1000

set.seed(0)
x <- rnorm(n, mean = -1, sd = 2)
mvn(modelName = "X", x) 

mu <- c(-1, 0, 1)

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% (2*diag(3)), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XII", x) 
mvn(modelName = "Spherical", x) 

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% diag(1:3), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XXI", x)
mvn(modelName = "Diagonal", x)

Sigma <- matrix(c(9,-4,1,-4,9,4,1,4,9), 3, 3)
set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) \%*\% chol(Sigma), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XXX", x) 
mvn(modelName = "Ellipsoidal", x) 
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
