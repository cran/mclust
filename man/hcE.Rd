\name{hcE}
\alias{hcE}
\alias{hcV}
\alias{hcEII}
\alias{hcVII}
\alias{hcEEE}
\alias{hcVVV}
\title{
  Model-based Hierarchical Clustering
}
\description{
  Agglomerative hierarchical clustering based on maximum likelihood
  for a Gaussian mixture model parameterized by eigenvalue decomposition.
}
\usage{
hcE(data, partition, minclus=1, \dots)
hcV(data, partition, minclus = 1, alpha = 1, \dots)
hcEII(data, partition, minclus = 1, \dots)
hcVII(data, partition, minclus = 1, alpha = 1, \dots)
hcEEE(data, partition, minclus = 1, \dots)
hcVVV(data, partition, minclus = 1, alpha = 1, beta = 1, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{partition}{
    A numeric or character vector representing a partition of
    observations (rows) of \code{data}. If provided, group merges will
    start with this partition. Otherwise, each observation is assumed to
    be in a cluster by itself at the start of agglomeration. 
  }
  \item{minclus}{
    A number indicating the number of clusters at which to stop the
    agglomeration. The default is to stop when all observations have been
    merged into a single cluster.
  }
  \item{alpha, beta}{
    Additional tuning parameters needed for initializatiion in some models. 
    For details, see Fraley 1998. The defaults provided are usually adequate.
  }
  \item{\dots}{
    Catch unused arguments from a \code{do.call} call.
  }
}
\value{
  A numeric two-column matrix in which the \emph{i}th row gives the minimum 
  index for observations in each of the two clusters merged at the
  \emph{i}th stage of agglomerative hierarchical clustering.
}
\details{
  Most models have memory usage of the order of the square of the
  number groups in the initial partition for fast execution.
  Some models, such as equal variance or \code{"EEE"},
  do not admit a fast algorithm under the usual agglomerative
  hierachical clustering paradigm. 
  These use less memory but are much slower to execute.
}
\section{References}{
  J. D. Banfield and A. E. Raftery (1993).
  Model-based Gaussian and non-Gaussian Clustering.
  \emph{Biometrics 49:803-821}. 
  
  C. Fraley (1998).
  Algorithms for model-based Gaussian hierarchical clustering.
  \emph{SIAM Journal on Scientific Computing 20:270-281}. 
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{hc}},
  \code{\link{hclass}}
}
\examples{
hcTree <- hcEII(data = iris[,-5])
cl <- hclass(hcTree,c(2,3))

\dontrun{
par(pty = "s", mfrow = c(1,1))
clPairs(iris[,-5],cl=cl[,"2"])
clPairs(iris[,-5],cl=cl[,"3"])

par(mfrow = c(1,2))
dimens <- c(1,2)
coordProj(iris[,-5], classification=cl[,"2"], dimens=dimens)
coordProj(iris[,-5], classification=cl[,"3"], dimens=dimens)
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
