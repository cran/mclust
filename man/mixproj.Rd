\name{mixproj}
\alias{mixproj}
\title{
Displays one standard deviation of an MVN mixture classification.
}
\section{DESCRIPTION}{
Displays one standard deviation of an MVN mxiture classification along
with data in for selected pairs of coordinates.
}
\usage{
mixproj(data, ms, partition, scale = F, newframe = T, k = 15, \dots)


}
\arguments{
\item{data}{
a matrix of onservations.
}
\item{ms}{
The result of an \code{mstep} calculation (a list consisting of \code{mu} and \code{sigma}
or \code{sigmasq}).
}
\item{partition}{
A integer vector giving an initail classification for each observation.
}
\item{dimens}{
A vector of length two giving the two variables of the data to be plotted.
}
\item{scale}{
A logical variable telling whether or not the same scale should be used for
both variables so as to preserve geometry. The default does not use the same 
scale.
}
\item{newframe}{
A logical variable indicating whether or not \code{frame} should be invoked before
plotting. The default is to invoke \code{frame}.
}
\item{k}{
Number of subdivisions for plotting segments of ellipsoids. Default: 8.
}
\item{...}{
use for the argument \code{symbols} indicating the desired symbols to be plotted 
(in the order that they appear in the classification).
}}
\seealso{
\code{\link{mstep}}, \code{\link{clpairs}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid = "VVV"),3)
z <- me( iris[,1:4], modelid = "VVV", ctoz(cl))
pars <- mstep(iris[,1:4], modelid="VVV", z)
mixproj(iris[,1:4], ms=pars, partition=ztoc(z), dimens=c(1,2))

}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.




