\name{mvn}
\alias{mvn}
\title{
  Univariate or Multivariate Normal Fit
}
\description{
  Computes the mean, covariance, and log-likelihood from fitting a single
  Gaussian to given data (univariate or multivariate normal).
}
\usage{
mvn( modelName, data, prior = NULL, warn = NULL, \dots)
}
\arguments{
  \item{modelName}{
    A character string representing a model name. This can be either
    \code{"Spherical"}, \code{"Diagonal"}, or \code{"Ellipsoidal"} or 
    else \cr
    \code{"X"} for one-dimensional data,\cr
    \code{"XII"} for a spherical Gaussian, \cr
    \code{"XXI"} for a diagonal Gaussian \cr
    \code{"XXX"} for a general ellipsoidal Gaussian 
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{prior}{
      Specification of a conjugate prior on the means and variances.
      The default assumes no prior.
  }
  \item{warn}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is given by \code{mclust.options("warn")}.
  }
 \item{\dots }{
   Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
 A list including the following components:
  \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{parameters}{
     \describe{
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}}
              for details.
        }
      }
  }
  \item{loglik}{
    The log likelihood for the data in the mixture model.
  }
  \item{Attributes:}{
      \code{"WARNING"} An appropriate warning if problems are 
      encountered in the computations.
  }
}
\seealso{
  \code{\link{mvnX}},
  \code{\link{mvnXII}},
  \code{\link{mvnXXI}},
  \code{\link{mvnXXX}},
  \code{\link{mclustModelNames}}
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
