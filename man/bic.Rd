\name{bic}
\alias{bic}
\title{
  BIC for Parameterized MVN Mixture Models
}
\description{
  Compute the BIC (Bayesian Information Criterion) for parameterized
  mixture models  given the loglikelihood, the dimension of the data,
  and number of mixture components in the model. 
}
\usage{
bic(modelName, loglik, n, d, G, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. Possible models: \cr \cr
    "E" for spherical, equal variance  (one-dimensional)  \cr
    "V" for spherical, variable variance (one-dimensional)  \cr
    "EII": spherical, equal volume  \cr
    "VII": spherical, unequal volume  \cr
    "EEI": diagonal, equal volume, equal shape  \cr
    "VEI": diagonal, varying volume, equal shape  \cr
    "EVI": diagonal, equal volume, varying shape  \cr
    "VVI": diagonal, varying volume, varying shape  \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation  \cr
    "EEV": ellipsoidal, equal volume and equal shape \cr
    "VEV": ellipsoidal, equal shape  \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
  }
  \item{loglik}{
    The loglikelihood for a data set with respect to the MVN mixture model
    specified in the \code{modelName} argument.
  }
  \item{n}{
    The number of observations in the data use to compute
    \code{loglik}.
  }
  \item{d}{
    The dimension of the data used to compute \code{loglik}.
  }
  \item{G}{
    The number of components in the MVN mixture model used to compute
    \code{loglik}.
  }
  \item{\dots}{
    Arguments for diagonal-specific methods, in particular
    \describe{
      \item{equalPro}{
	A logical variable indicating whether or not the components in the
	model are assumed to be present in equal proportion. The default is
	\code{.Mclust\$equalPro}. 
      }
      \item{noise}{
	A logical variable indicating whether or not the model includes and
	optional Poisson noise component. The default is to assume that the
	model does not include a noise component.
      }
    }
  }
}
\value{
  The BIC or Bayesian Information Criterion for the given input arguments.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{bicE}}, \ldots,
  \code{\link{bicVVV}},
  \code{\link{EMclust}},
  \code{\link{estep}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}.
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

n <- nrow(irisMatrix)
d <- ncol(irisMatrix)
G <- 3

emEst <- me(modelName="VVI",data=irisMatrix, unmap(irisClass))
names(emEst)

args(bic)
bic(modelName="VVI",loglik=emEst$loglik,n=n,d=d,G=G)
do.call("bic", emEst)    ## alternative call
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
