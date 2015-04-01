\name{nVarParams}
\alias{nVarParams}
\title{
  Number of Variance Parameters in Gaussian Mixture Models
}
\description{
 Gives the number of variance parameters for parameterizations of the 
 Gaussian mixture model that are used in MCLUST. 
}
\usage{
nVarParams(modelName, d, G, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
  \item{d}{
    The dimension of the data. Not used for models in which neither
    the shape nor the orientation varies.
  }
  \item{G}{
    The number of components in the Gaussian mixture model used to compute
    \code{loglik}.
  }
  \item{\dots }{
   Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  The number of variance parameters in the corresponding Gaussian mixture
  model. 
}
\details{
  To get the total number of parameters in model, add \code{G*d} for the
  means and \code{G-1} for the mixing proportions if they are unequal.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
 
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{bic}}, \code{\link{nMclustParams}}.
}
\examples{
mapply(nVarParams, mclust.options("emModelNames"), d = 2, G = 3)
}
\keyword{cluster}
