\name{Mclust}
\alias{Mclust}
\title{
  Model-Based Clustering
}
\description{
  Clustering via EM initialized by hierarchical clustering for
  parameterized Gaussian mixture models. The number of clusters and the
  clustering model is chosen to maximize the BIC. 
}
\usage{
Mclust(data, minG, maxG)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{minG}{
    An integer vector specifying the minimum number of mixture components 
    (clusters) to be considered.
    The default is \code{1} component.
  }
  \item{maxG}{
    An integer vector specifying the maximum number of mixture components 
    (clusters) to be considered.
    The default is \code{9} components.
  }
}
\value{
  A list representing the best model (according to BIC) for the given range
  of numbers of clusters. The following components are included: 
  \item{BIC}{ 
    A matrix giving the BIC value for each model (rows) and number of
    clusters (columns).
  }
  \item{bic}{ 
    A scalar giving the optimal BIC value.
  }
  \item{modelName}{ 
    The MCLUST name for the best model according to BIC. 
  }
  \item{classification}{ 
    The classification corresponding to the optimal BIC value. 
  }
  \item{uncertainty}{ 
    The  uncertainty in the classification corresponding to the
    optimal BIC value.
  }
  \item{mu}{ 
    For multidimensional models, a matrix whose columns are the means
    of each group in the best model. For one-dimensional models, a
    vector whose entries are the means for each group in the best
    model. 
  }
  \item{sigma}{ 
    For multidimensional models, a three dimensional array in which
    \code{sigma[,,k]} gives the covariance for the \emph{k}th group in
    the best model. For one-dimensional models, either a scalar giving
    a common variance for the groups or a vector whose entries are the
    variances for each group in the best model.
  }
  \item{pro}{ 
    The mixing probabilities for each component in the best model.
  }
  \item{z}{ 
    A matrix whose \emph{[i,k]}th entry is the probability that
    observation \emph{i} belongs to the \emph{k} component in the
    model. The optimal classification is derived from this, chosing
    the class to be the one giving the maximum probability.
  }
  \item{loglik}{ 
    The log likelihood for the data under the best model.
  }
}
\section{Details}{
  The following models are compared in \code{Mclust}: \cr\cr 
  "E" for spherical, equal variance (one-dimensional) \cr
  "V" for spherical, variable variance (one-dimensional) \cr\cr
  "EII": spherical, equal volume \cr
  "VII": spherical, unequal volume \cr
  "EEI": diagonal, equal volume, equal shape \cr
  "VVI": diagonal, varying volume, varying shape \cr
  "EEE": ellipsoidal, equal volume, shape, and orientation \cr
  "VVV": ellipsoidal, varying volume, shape, and orientation \cr\cr
  \code{Mclust} is intended to combine \code{EMclust} and its
  \code{summary} in a simiplified one-step model-based clustering
  function. The latter provide more flexibility including choice of
  models. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}.
 
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{plot.Mclust}},
  \code{\link{EMclust}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]
irisMclust <- Mclust(irisMatrix)

\dontrun{plot(irisMclust,irisMatrix)}
}
\keyword{cluster}
% docclass is function
