\name{bic}
\alias{bic}
\title{
  BIC for Parameterized Gaussian Mixture Models
}
\description{
  Computes the BIC (Bayesian Information Criterion) for parameterized
  mixture models given the loglikelihood, the dimension of the data,
  and number of mixture components in the model. 
}
\usage{
bic(modelName, loglik, n, d, G, noise=FALSE, equalPro=FALSE, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
  \item{loglik}{
    The log-likelihood for a data set with respect to the Gaussian mixture model
    specified in the \code{modelName} argument.
  }
  \item{n}{
    The number of observations in the data used to compute \code{loglik}.
  }
  \item{d}{
    The dimension of the data used to compute \code{loglik}.
  }
  \item{G}{
    The number of components in the Gaussian mixture model used to compute
    \code{loglik}.
  }
 \item{noise}{
   A logical variable indicating whether or not the model includes an
   optional Poisson noise component. The default is to assume no noise
   component.
  }
 \item{equalPro}{
    A logical variable indicating whether or not the components in the
    model are assumed to be present in equal proportion. The default is
    to assume unequal mixing proportions.
  }
  \item{\dots}{
    Catches unused arguments in an indirect or list call via \code{do.call}.
  }
}
\value{
  The BIC or Bayesian Information Criterion for the given input arguments.
}

\seealso{
  \code{\link{mclustBIC}},
  \code{\link{nVarParams}},
  \code{\link{mclustModelNames}}.
}
\examples{
\dontrun{
n <- nrow(iris)
d <- ncol(iris)-1
G <- 3

emEst <- me(modelName="VVI", data=iris[,-5], unmap(iris[,5]))
names(emEst)

args(bic)
bic(modelName="VVI", loglik=emEst$loglik, n=n, d=d, G=G)
# do.call("bic", emEst)    ## alternative call
}
}
\keyword{cluster}
