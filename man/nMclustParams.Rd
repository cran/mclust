\name{nMclustParams}
\alias{nMclustParams}

\title{Number of Estimated Parameters in Gaussian Mixture Models}

\description{
 Gives the number of estimated parameters for parameterizations of the 
 Gaussian mixture model that are used in MCLUST. 
}
\usage{
nMclustParams(modelName, d, G, noise = FALSE, equalPro = FALSE, \dots)
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
  \item{noise}{
    A logical variable indicating whether or not the model includes an
   optional Poisson noise component. 
  }
  \item{equalPro}{
    A logical variable indicating whether or not the components in the
    model are assumed to be present in equal proportion.
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
\seealso{
  \code{\link{bic}}, \code{\link{nVarParams}}.
}
\examples{
mapply(nMclustParams, mclust.options("emModelNames"), d = 2, G = 3)
}
\keyword{cluster}
