\name{mvnX}
\alias{mvnX}
\alias{mvnXII}
\alias{mvnXXI}
\alias{mvnXXX}
\title{
  Univariate or Multivariate Normal Fit
}
\description{
  Computes the mean, covariance, and log-likelihood from fitting a single
  Gaussian (univariate or multivariate normal).
}
\usage{
mvnX(data, prior = NULL, warn = NULL, \dots)
mvnXII(data, prior = NULL, warn = NULL, \dots)
mvnXXI(data, prior = NULL, warn = NULL, \dots)
mvnXXX(data, prior = NULL, warn = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
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

\details{
\describe{
\item{\code{mvnXII}}{computes the best fitting Gaussian with the covariance restricted to be a multiple of the identity.}
\item{\code{mvnXXI}}{computes the best fitting Gaussian with the covariance restricted to be diagonal.}
\item{\code{mvnXXX}}{computes the best fitting Gaussian with ellipsoidal (unrestricted) covariance.}
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
  \code{\link{mvn}},
  \code{\link{mstepE}}
}
\examples{
\donttest{
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
}
\keyword{cluster}
