\name{mhclass}
\alias{mhclass}
\title{
Classifications for hierarchical clustering.
}
\description{
Gives classifications corresponding to the stages of hierarchical clustering.
}
\usage{
mhclass(tree, nclusters)


}
\arguments{
\item{tree}{
An \code{"mhtree"} object, which is a matrix whose columns are integers 
corresponding to clusters merged at each stage of hierarchical clustering
(classification tree), together with an attribute giving the initial
partition for the hierarchical process.
}
\item{nclusters}{
the numbers of clusters for which classifications are desired. The default
is to produce all classifications implied by \code{tree} given the initial
partition.
}}
\value{
a matrix whose columns are the classifications corresponding to each value of
\code{nclusters}. Columns are labeled by the \code{as.character} representation of
the number of clusters.
}
\seealso{
\code{\link{mhtree}}
}
\examples{
data(iris)
mhclass(mhtree(iris[,1:4]), 2:6)
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.




