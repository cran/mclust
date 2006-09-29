\name{priorControl}
\alias{priorControl}
\title{
  Conjugate Prior for Gaussian Mixtures.
}
\description{
   Specify a conjugate prior for Gaussian mixtures.
}
\usage{
priorControl(functionName = "defaultPrior", \dots) 
}
\arguments{
  \item{functionName}{
    The name of the function specifying the conjugate prior.
    The default function is \code{defaultPrior}, which can be
    used a template for alternative specification.  
  }
  \item{\dots}{
     Optional named arguments to the function specified in \code{functionName}
     together with their values.
    }
  }
}
\value{
  A list with the function name as the first component. The remaining
  components (if any) consist of a list of arguments to the function
  with assigned values.
}
\details{
  \code{priorControl} is used to specify a conjugate prior  
  for EM within \emph{MCLUST}.
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
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{mclustBIC}},
  \code{\link{me}},
  \code{\link{mstep}},
  \code{\link{defaultPrior}}
}
\examples{
# default prior
irisBIC <- mclustBIC(iris[,-5], prior = priorControl())
summary(irisBIC, iris[,-5])

# no prior on the mean; default prior on variance
irisBIC <- mclustBIC(iris[,-5], prior = priorControl(shrinkage = 0))
summary(irisBIC, iris[,-5])
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
