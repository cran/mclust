\name{mixproj}
\alias{mixproj}
\alias{mvn2plot}
\title{
Displays one standard deviation of an MVN mixture classification.
}
\description{
\code{mixproj} displays one standard deviation of an MVN mxiture
classification along with data in for selected pairs of
coordinates. Function \code{mvn2plot} is used to plot the confidence
ellipses.
}
\usage{
mixproj(data, ms, partition, dimens, scale = F, k = 15, xlim, ylim,
xlab, ylab, col=partition, pch=partition, \dots)
mvn2plot(mu, sigma, k, add=T, col=NULL, xlim=NULL, ylim=NULL, xlab=NULL,
ylab=NULL, \dots)
}

\arguments{
\item{data}{
a matrix of onservations.
}
\item{ms}{
The result of an \code{mstep} calculation (a list consisting of
\code{mu} and \code{sigma}).
}
\item{partition}{
A integer vector giving an initial classification for each observation.
}
\item{dimens}{
  A vector of length two giving the two variables of the data to be
  plotted. If not provided, a random selection of two variables is made.
  
}
\item{scale}{
A logical variable telling whether or not the same scale should be used for
both variables so as to preserve geometry. The default does not use the same 
scale.
}
\item{k}{
Number of subdivisions for plotting segments of ellipsoids. Default: 8.
}
\item{xlim, ylim, xlab, ylab, col, pch}{
Standard arguments to \code{plot} with default values.
}
\item{add}{Add ellipses from \code{mvn2plot} to the existing plot.}

\item{...}{
  Further arguments to \code{plot}.
}
\item{mu}{
  Cluster centers.
}
\item{sigma}{
  Cluster covariance matrices.}}
\seealso{
\code{\link{mstep}}, \code{\link{clpairs}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid = "VVV"),3)
z <- me( iris[,1:4], modelid = "VVV", ctoz(cl))
pars <- mstep(iris[,1:4], modelid="VVV", z)
mixproj(iris[,1:4], ms=pars, partition=ztoc(z), dimens=c(1,2))
### Plot only one cluster
mvn2plot(pars$mu[1:2,1], pars$sigma[1:2, 1:2, 1], 25, F, col="red")
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.




