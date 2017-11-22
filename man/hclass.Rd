\name{hclass}
\alias{hclass}
\title{
  Classifications from Hierarchical Agglomeration 
}
\description{
  Determines the classifications corresponding to different numbers of groups
  given merge pairs from hierarchical agglomeration.
}
\usage{
hclass(hcPairs, G)
}
\arguments{
  \item{hcPairs}{
    A numeric two-column matrix in which the \emph{i}th row gives the minimum 
    index for observations in each of the two clusters merged at the
    \emph{i}th stage of agglomerative hierarchical clustering.
  }
  \item{G}{
    An integer or vector of integers giving the number of clusters for which
    the corresponding classfications are wanted.
  }
}
\value{
  A matrix with \code{length(G)} columns, each column 
  corresponding to a classification. Columns are indexed by the character
  representation of the integers in \code{G}.
}
\seealso{
  \code{\link{hc}},
  \code{\link{hcE}}
}
\examples{
hcTree <- hc(modelName="VVV", data = iris[,-5])
cl <- hclass(hcTree,c(2,3))

\dontrun{
par(pty = "s", mfrow = c(1,1))
clPairs(iris[,-5],cl=cl[,"2"])
clPairs(iris[,-5],cl=cl[,"3"])
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
