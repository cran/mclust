\name{sigma2decomp}
\alias{sigma2decomp}
\title{
  Convert mixture component covariances to decomposition form.
}
\description{
  Converts a set of covariance matrices from representation as a 3-D array 
  to a parameterization by eigenvalue decomposition.
}
\usage{
sigma2decomp(sigma, G, tol, \dots)
}
\arguments{
  \item{sigma}{
    Either a 3-D array whose [,,k]th component is the covariance matrix for the
    kth component in an MVN mixture model, or a single covariance
    matrix in the case that all components have the same covariance.
  }
  \item{G}{
    The number of components in the mixture. When 
    \code{sigma} is a 3-D array, the number of components
    can be inferred from its dimensions.
  }
  \item{tol}{
    Tolerance for determining whether or not the covariances have equal volume,
    shape, and or orientation. The default is the square root of the relative
    machine precision, \code{sqrt(.Machine\$double.eps)}, 
    which is about \code{1.e-8}.
  }
  \item{\dots}{
    Catch unused arguments from a \code{do.call} call.
  }
}
\value{
  The covariance matrices for the mixture components in decomposition form,
  including the following components: 
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
    Either a \emph{G} by \emph{d} matrix in which the \emph{k}th
    column is the shape of the covariance matrix (normalized to have
    determinant 1) for the \emph{k}th component, or a \emph{d}-vector
    giving a common shape for all components. 
  }
  \item{orientation}{
    Either a \emph{d} by \emph{d} by \emph{G} array whose
    \code{[,,k]}th entry is the orthonomal matrix of eigenvectors of
    the covariance matrix of the \emph{k}th component, or a \emph{d}
    by \emph{d} orthonormal matrix if the mixture components have a
    common orientation. The \code{orientation} component of
    \code{decomp} can be omitted in spherical and diagonal models, for
    which the principal components are parallel to the coordinate axes
    so that the orientation matrix is the identity.  
  }
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation, and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{decomp2sigma}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

meEst <- meEEE(irisMatrix, unmap(irisClass)) 
names(meEst)
meEst$sigma

sigma2decomp(meEst$sigma)
do.call("sigma2decomp", meEst)  ## alternative call
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
