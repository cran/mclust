\name{clpairs}
\alias{clpairs}
\title{
Classifications for hierarchical clustering.
}
\section{DESCRIPTION}{
Gives classifications corresponding to the stages of hierarchical clustering.
}
\usage{
clpairs(x, partition, symbols, two=T, newframe = T, \dots)


}
\arguments{
\item{x}{
matrix of observations.
}
\item{partition}{
a vector of integers indicating a partition of the data.
}
\item{symbols}{
graphics symbols denoting the different groups in the partition (in order
of appearance).
}
\item{two}{
logical variable indicating whether or not the two group case should be
displayed as a single plot or a pairs plot. The default \code{two = T} displays 
is as a single plot.
}
\item{newframe}{
logical variable indicating whether or not \code{frame} should be invoked before
plotting. The default is to invoke \code{frame}.	
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
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.




