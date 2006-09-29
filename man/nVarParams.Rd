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
nVarParams(modelName, d, G)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{mclustModelNames} describes the available models.
  }
  \item{d}{
    The dimension of the data. Not used for models in which neither
    the shape nor the orientation varies.
  }
  \item{G}{
    The number of components in the Gaussian mixture model used to compute
    \code{loglik}.
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
 
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{bic}}
}
\examples{
sapply(.Mclust$emModelNames, nVarParams, d=2, G=1)
}
\keyword{cluster}
