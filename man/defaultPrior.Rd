\name{defaultPrior}
\alias{defaultPrior}
\title{
  Default conjugate prior for Gaussian mixtures.
}
\description{
   Default conjugate prior specification for Gaussian mixtures.
}
\usage{
defaultPrior(data, G, modelName, \dots) 
}
\arguments{
  \item{data}{
    The name of the function specifying the conjgate prior.
    The default function is \code{defaultPrior}, which can be
    used a template for  
  }
  \item{G}{
    The number of mixture components.
  }
  \item{modelName}{
    A character string indicating the model: \cr
    "E": equal variance  (one-dimensional) \cr
    "V": variable variance (one-dimensional)\cr 
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape\cr 
    "VEI": diagonal, varying volume, equal shape\cr 
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
  }
  \item{\dots}{
     One or more of the following:
     \itemize{
      \item{dof}{
        The degrees of freedom for the prior on the variance. 
        The default is \code{d + 2}, where \code{d} is
        the dimension of the data.
      }
      \item{scale}{
        The scale parameter for the prior on the variance. 
        The default is \code{var(data)/G^(2/d)},
        where \code{d} is the domension of the data.
      }
      \item{shrinkage}{
        The shrinkage parameter for the prior on the mean. 
        The default value is 0.01. 
        If 0 or NA, no prior is assumed for the mean.
      }
      \item{mean}{
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
  \code{defaultPrior} is as a default prior specification  
  for EM within \emph{MCLUST}. It is usually not necessary
  to invoke \code{defaultPrior} explicitly (it does not appear in the examples
  below because it is the default function name in \code{priorControl}).
  This function allows considerable flexibility in the prior specification,
  and can be used as a template for further users that want to specify
  their own conjugate prior beyond what the arguments will allow. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2005).
  Bayesian regularization for normal mixture estimation and model-based
  clustering.
  Technical Report, Department of Statistics, University of Washington.

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
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
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
