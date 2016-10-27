\name{RARI}
\alias{RARI}
\title{
  Ranked Adjusted Rand Index
}
\description{
   Computes the ranked adjusted Rand index comparing two classifications with supplementary distance matrices. 
}
\usage{
RARI(x, y, dist_x, dist_y)
}
\arguments{
  \item{x}{
    A numeric or character vector of class labels.
  }
  \item{y}{
    A numeric or character vector of class labels.
    The length of \code{y} should be the same as that of \code{x}.
  }
  \item{dist_x}{
    A distance matrix of class 'dist' or 'matrix' for \code{x} class labels.
  }
  \item{dist_y}{
    A distance matrix of class 'dist' or 'matrix' for \code{y} class labels.
  }
}
\value{
  The ranked adjusted Rand index comparing the two partitions and their distance matrices (a scalar).  
  This index has zero expected value only in the case of (a) random partition, (b) when in one clustering
  all entities are in the same clustering and in the other, every entitity is in its own cluster, and (c)
  all clusters are equidistant from each other. The ranked adjusted Rand index is bounded above by 1 in the
  case of (a) perfect agreement between two partitions and (b) equally ranked relative distances between clusters.
  
  If the arguments \code{dist_x} and \code{dist_y} are \code{NULL}, the output is identical to the adjusted Rand index.
}

\section{References}{
 Pinto, F.R., Carrico, J.A., Ramirez, M., and Almeida, J.S. (2007). Ranked Adjusted
     Rand: integrating distance and partition information in a measure of clustering
     agreement. \emph{BMC Bioinformatics 8:44}. doi: 10.1186/1471-2105-8-44
}
\seealso{
  \code{\link{classError}},
  \code{\link{mapClass}},
  \code{\link{table}}
}
\examples{
data(iris)
x <- iris$Species
y <- iris$Species
dist_x <- dist(iris[, 1:4]) # Using all measures in the distance matrix
dist_y <- dist(iris[, 2])   # Using only Sepal.Width in the distance matrix
  
adjustedRandIndex(x, y)
RARI(x, y, dist_x, dist_y)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
