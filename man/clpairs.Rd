\name{clpairs}
\alias{clpairs}
\title{
Classifications for hierarchical clustering.
}
\description{
Gives classifications corresponding to the stages of hierarchical clustering.
}
\usage{
clpairs(x, partition, col=partition, \dots)
}

\arguments{
\item{x}{
matrix of observations.
}
\item{partition}{
a vector of integers indicating a partition of the data.
}
\item{col}{
the default is to color the objects according to the partition.
}
\item{...}{
other arguments as specified by \code{pairs}
}}
\value{
a matrix whose columns are the classifications corresponding to each value of
\code{nclusters}.
}
\seealso{
\code{\link{pairs}}, \code{\link{mixproj}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid="VVV"), 3)
clpairs(iris[,1:4],cl)
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.




