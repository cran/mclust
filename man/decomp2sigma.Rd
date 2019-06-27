\name{decomp2sigma}
\alias{decomp2sigma}
\title{
  Convert mixture component covariances to matrix form
}
\description{
  Converts covariances from a parameterization by  eigenvalue decomposition 
  or cholesky factorization to representation as a 3-D array.
}
\usage{
decomp2sigma(d, G, scale, shape, orientation, \dots)
}
\arguments{
  \item{d}{
    The dimension of the data. 
  }
  \item{G}{
    The number of components in the mixture model. 
  }
  \item{scale}{
    Either a \emph{G}-vector giving the scale of the covariance (the
    \emph{d}th root of its determinant) for each component in the
    mixture model, or a single numeric value if the scale is the same
    for each component.  
  }
  \item{shape}{
    Either a \emph{G} by \emph{d} matrix in which the \emph{k}th column
    is the shape of the covariance matrix (normalized to have
    determinant 1) for the \emph{k}th component, or a \emph{d}-vector
    giving a common shape for all components. 
  }
  \item{orientation}{
    Either a \emph{d} by \emph{d} by \emph{G} array whose \code{[,,k]}th
    entry is the orthonomal matrix whose columns are the eigenvectors
    of the covariance matrix of the \emph{k}th component, or a
    \emph{d} by \emph{d} orthonormal matrix if the mixture components have a common
    orientation. The \code{orientation} component of \code{decomp} can
    be omitted in spherical and diagonal models, for which the principal
    components are parallel to the coordinate axes so that the
    orientation matrix is the identity.  
  }
  \item{\dots}{
    Catches unused arguments from an indirect or list call via \code{do.call}.
  }
}
\value{
  A 3-D array whose \code{[,,k]}th component is the 
  covariance matrix of the \emph{k}th component in an MVN mixture model.
}

\seealso{
  \code{\link{sigma2decomp}}
}
\examples{
meEst <- meVEV(iris[,-5], unmap(iris[,5])) 
names(meEst)
meEst$parameters$variance

dec <- meEst$parameters$variance
decomp2sigma(d=dec$d, G=dec$G, shape=dec$shape, scale=dec$scale,
             orientation = dec$orientation)
\dontrun{
do.call("decomp2sigma", dec)  ## alternative call
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
