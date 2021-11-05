\name{hc}
\alias{hc}
\alias{print.hc}

\title{Model-based Agglomerative Hierarchical Clustering}

\description{
  Agglomerative hierarchical clustering based on maximum likelihood criteria 
  for Gaussian mixture models parameterized by eigenvalue decomposition.
}

\usage{
hc(data,
   modelName = "VVV",  
   use = "VARS",
   partition = dupPartition(data), 
   minclus = 1, \dots)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations (\eqn{n}) and
    columns correspond to variables (\eqn{d}).
  }
  \item{modelName}{
    A character string indicating the model to be used in model-based agglomerative hierarchical clustering.\cr
    Possible models are:
    \describe{
    \item{\code{"E"}}{equal variance (one-dimensional);}
    \item{\code{"V"}}{spherical, variable variance (one-dimensional);}
    \item{\code{"EII"}}{spherical, equal volume;}
    \item{\code{"VII"}}{spherical, unequal volume;}
    \item{\code{"EEE"}}{ellipsoidal, equal volume, shape, and orientation;}
    \item{\code{"VVV"}}{ellipsoidal, varying volume, shape, and orientation (default).}
    }
    If \code{hc()} is used for initialization of EM algorithm then the default is taken from \code{mclust.options("hcModelName")}. See \code{\link{mclust.options}}.
  }
  \item{use}{
    A character string specifying the type of input variables/data transformation to be used for model-based agglomerative hierarchical clustering.\cr
    Possible values are:
    \describe{
    \item{\code{"VARS"}}{original variables (default);}
    \item{\code{"STD"}}{standardized variables (centered and scaled);}
    \item{\code{"SPH"}}{sphered variables (centered, scaled and uncorrelated)  
    computed using SVD;}
    \item{\code{"PCS"}}{principal components computed using SVD on centered 
    variables (i.e. using the covariance matrix);}
    \item{\code{"PCR"}}{principal components computed using SVD on standardized 
    (center and scaled) variables (i.e. using the correlation matrix);}
    \item{\code{"SVD"}}{scaled SVD transformation.}
    }
    If \code{hc()} is used for initialization of EM algorithm then the default is taken from \code{mclust.options("hcUse")}. See \code{\link{mclust.options}}.\cr
    For further details see Scrucca and Raftery (2015).
  }
  \item{partition}{
    A numeric or character vector representing a partition of
    observations (rows) of \code{data}. 
    If provided, group merges will start with this partition. 
    Otherwise, each observation is assumed to be in a cluster by itself 
    at the start of agglomeration.
    Starting with version 5.4.8, by default the function
    \code{\link{dupPartition}} is used to start with all duplicated
    observations in the same group, thereby keeping duplicates in the 
    same group throughout the modelling process.
  }
  \item{minclus}{
    A number indicating the number of clusters at which to stop the
    agglomeration. The default is to stop when all observations have been
    merged into a single cluster.
  }
  \item{\dots}{
    Arguments for the method-specific \code{hc} functions. See for example
    \code{\link{hcE}}.
  }
}

\value{
  The function \code{hc()} returns a numeric two-column matrix in which 
  the \emph{i}th row gives the minimum index for observations in each of 
  the two clusters merged at the \emph{i}th stage of agglomerative 
  hierarchical clustering. Several other informations are also returned
  as attributes.
  
  The plotting method \code{plot.hc()} draws a dendrogram, which can be based
  on either the classification loglikelihood or the merge level (number of
  clusters). For details, see the associated help file.
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
  covariances), then underlying model is the same as that for
  Ward's method for hierarchical clustering.
}
\references{
  Banfield J. D. and Raftery A. E. (1993).
  Model-based Gaussian and non-Gaussian Clustering.
  \emph{Biometrics}, 49:803-821. 
  
  Fraley C. (1998).
  Algorithms for model-based Gaussian hierarchical clustering.
  \emph{SIAM Journal on Scientific Computing}, 20:270-281. 
  
  Fraley C. and Raftery A. E. (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association}, 97:611-631. 
  
  Scrucca L. and Raftery A. E. (2015).
  Improved initialisation of model-based clustering using Gaussian hierarchical partitions. 
  \emph{Advances in Data Analysis and Classification}, 9/4:447-460.
}

\seealso{
  \code{\link{hcE}}, \dots,
  \code{\link{hcVVV}},
  \code{\link{plot.hc}},
  \code{\link{hclass}},
  \code{\link{mclust.options}}  
}

\examples{
hcTree <- hc(modelName = "VVV", data = iris[,-5])
hcTree
cl <- hclass(hcTree,c(2,3))
table(cl[,"2"])
table(cl[,"3"])

\donttest{
clPairs(iris[,-5], classification = cl[,"2"])
clPairs(iris[,-5], classification = cl[,"3"])
}
}
\keyword{cluster}
