\name{hc}
\alias{hc}
\alias{print.hc}
\alias{plot.hc}
\alias{as.dendrogram.hc}
\alias{as.hclust.hc}

\title{Model-based Agglomerative Hierarchical Clustering}

\description{
  Agglomerative hierarchical clustering based on maximum likelihood criteria 
  for Gaussian mixture models parameterized by eigenvalue decomposition.
}

\usage{
hc(data,
   modelName = mclust.options("hcModelName"),  
   use = mclust.options("hcUse"), \dots)

\method{plot}{hc}(x, \dots)

\method{as.dendrogram}{hc}(object, \dots)

\method{as.hclust}{hc}(x, \dots)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations (\eqn{n}) and
    columns correspond to variables (\eqn{d}).
  }
  \item{modelName}{
    A character string indicating the model to be used.\cr
    Possible models are:
    \describe{
    \item{\code{"E"}}{equal variance (one-dimensional)}
    \item{\code{"V"}}{spherical, variable variance (one-dimensional)}
    \item{\code{"EII"}}{spherical, equal volume}
    \item{\code{"VII"}}{spherical, unequal volume}
    \item{\code{"EEE"}}{ellipsoidal, equal volume, shape, and orientation}
    \item{\code{"VVV"}}{ellipsoidal, varying volume, shape, and orientation.}
    }
    By default the model provided by \code{mclust.options("hcModelName")} is used. 
    See \code{\link{mclust.options}}.
  }
  \item{use}{
    A string or a vector of character strings specifying the type of input 
    variables/data transformation to be used for model-based hierarchical 
    clustering.\cr
    By default the method specified in \code{mclust.options("hcUse")} is used. 
    See \code{\link{mclust.options}}.
  }
  \item{\dots}{
    Arguments for the method-specific \code{hc} functions. See for example
    \code{\link{hcE}}.
  }
  \item{object, x}{
    An object of class \code{'hc'} resulting from a call to \code{hc()}.
  }

}
\value{
  The function \code{hc()} returns a numeric two-column matrix in which 
  the \emph{i}th row gives the minimum index for observations in each of 
  the two clusters merged at the \emph{i}th stage of agglomerative 
  hierarchical clustering. Several other informations are also returned
  as attributes.
  
  The plotting function \code{plot.hc()} draws a dendrogram by first
  converting the input object from class \code{'hc'} to class 
  \code{'dendrogram'} and then plot the transformed object using 
  \code{\link{plot.dendrogram}}.
  
  The functions \code{as.dendrogram.hc()} and \code{as.hclust.hc()} are
  used to convert the input object from class \code{'hc'} to class, 
  respectively, \code{'dendrogram'} and \code{'hclust'}.
}

\details{
  Most models have memory usage of the order of the square of the
  number groups in the initial partition for fast execution.
  Some models, such as equal variance or \code{"EEE"},
  do not admit a fast algorithm under the usual agglomerative
  hierarchical clustering paradigm. 
  These use less memory but are much slower to execute.
}

\note{
  If \code{modelName = "E"} (univariate with equal variances) or
  \code{modelName = "EII"} (multivariate with equal spherical
  covariances), then the method is equivalent to Ward's method for
  hierarchical clustering.
}
\references{
  J. D. Banfield and A. E. Raftery (1993).
  Model-based Gaussian and non-Gaussian Clustering.
  \emph{Biometrics 49:803-821}. 
  
  C. Fraley (1998).
  Algorithms for model-based Gaussian hierarchical clustering.
  \emph{SIAM Journal on Scientific Computing 20:270-281}. 
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
}
\seealso{
  \code{\link{hcE}},...,
  \code{\link{hcVVV}},
  \code{\link{hclass}},
  \code{\link{mclust.options}}  
}
\examples{
hcTree <- hc(modelName = "VVV", data = iris[,-5])
cl <- hclass(hcTree,c(2,3))

\dontrun{
par(pty = "s", mfrow = c(1,1))
clPairs(iris[,-5],cl=cl[,"2"])
clPairs(iris[,-5],cl=cl[,"3"])

par(mfrow = c(1,2))
dimens <- c(1,2)
coordProj(iris[,-5], dimens = dimens, classification=cl[,"2"])
coordProj(iris[,-5], dimens = dimens, classification=cl[,"3"])
}
}
\keyword{cluster}
