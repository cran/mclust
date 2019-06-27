\name{defaultPrior}
\alias{defaultPrior}
\title{
  Default conjugate prior for Gaussian mixtures
}
\description{
   Default conjugate prior specification for Gaussian mixtures.
}
\usage{
defaultPrior(data, G, modelName, \dots) 
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{G}{
    The number of mixture components.
  }
  \item{modelName}{
    A character string indicating the model: \cr\cr
    \code{"E"}: equal variance  (univariate) \cr
    \code{"V"}: variable variance (univariate)\cr 
    \code{"EII"}: spherical, equal volume \cr
    \code{"VII"}: spherical, unequal volume \cr
    \code{"EEI"}: diagonal, equal volume and shape\cr 
    \code{"VEI"}: diagonal, varying volume, equal shape\cr 
    \code{"EVI"}: diagonal, equal volume, varying shape \cr
    \code{"VVI"}: diagonal, varying volume and shape \cr
    \code{"EEE"}: ellipsoidal, equal volume, shape, and orientation \cr
    \code{"EEV"}: ellipsoidal, equal volume and equal shape\cr
    \code{"VEV"}: ellipsoidal, equal shape \cr
    \code{"VVV"}: ellipsoidal, varying volume, shape, and orientation. \cr\cr
    A description of the models above is provided in the help of 
    \code{\link{mclustModelNames}}. Note that in the multivariate case 
    only 10 out of 14 models may be used in conjunction with a prior, i.e.
    those available in \emph{MCLUST} up to version 4.4.
  }
  \item{\dots}{
     One or more of the following:
     \describe{
      \item{\code{dof}}{
        The degrees of freedom for the prior on the variance. 
        The default is \code{d + 2}, where \code{d} is
        the dimension of the data.
      }
      \item{\code{scale}}{
        The scale parameter for the prior on the variance. 
        The default is \code{var(data)/G^(2/d)},
        where \code{d} is the dimension of the data.
      }
      \item{\code{shrinkage}}{
        The shrinkage parameter for the prior on the mean. 
        The default value is 0.01. 
        If 0 or NA, no prior is assumed for the mean.
      }
      \item{\code{mean}}{
        The mean parameter for the prior. 
        The default value is \code{colMeans(data)}.
      }                   
    }
  }
}
\value{
  A list giving the prior degrees of freedom, scale, shrinkage, and mean.
}
\details{
  \code{defaultPrior} is a function whose default is to output the
  default prior specification for EM within \emph{MCLUST}.\cr
  Furthermore, \code{defaultPrior} can be used as a template to specify 
  alternative parameters for a conjugate prior.
}
\references{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association} 97:611-631. 

  C. Fraley and A. E. Raftery (2005, revised 2009).
  Bayesian regularization for normal mixture estimation and model-based
  clustering.
  Technical Report, Department of Statistics, University of Washington.

  C. Fraley and A. E. Raftery (2007).
  Bayesian regularization for normal mixture estimation and model-based
  clustering. \emph{Journal of Classification} 24:155-181.
}
\seealso{
  \code{\link{mclustBIC}},
  \code{\link{me}},
  \code{\link{mstep}},
  \code{\link{priorControl}}
}
\examples{
# default prior
irisBIC <- mclustBIC(iris[,-5], prior = priorControl())
summary(irisBIC, iris[,-5])

# equivalent to previous example
irisBIC <- mclustBIC(iris[,-5], 
                     prior = priorControl(functionName = "defaultPrior"))
summary(irisBIC, iris[,-5])

# no prior on the mean; default prior on variance
irisBIC <- mclustBIC(iris[,-5], prior = priorControl(shrinkage = 0))
summary(irisBIC, iris[,-5])

# equivalent to previous example
irisBIC <- mclustBIC(iris[,-5], prior =
                     priorControl(functionName="defaultPrior", shrinkage=0))
summary(irisBIC, iris[,-5])

defaultPrior( iris[-5], G = 3, modelName = "VVV")
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
