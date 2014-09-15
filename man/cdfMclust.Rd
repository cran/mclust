\name{cdfMclust}
\alias{cdfMclust}
\alias{quantileMclust}

\title{
  Cumulative Distribution and Quantiles for a univariate Gaussian mixture 
  distribution
}

\description{
Compute the cumulative density function (cdf) or quantiles from an estimated one-dimensional Gaussian mixture fitted using \code{\link{densityMclust}}.}

\usage{
cdfMclust(object, data, ngrid = 100, \dots)
quantileMclust(object, p, ...)
}

\arguments{
  \item{object}{a \code{densityMclust} model object.}
  \item{data}{a numeric vector of evaluation points.}
  \item{ngrid}{the number of points in a regular grid to be used as evaluation points if no \code{data} are provided.}
  \item{p}{a numeric vector of probabilities.}
 \item{\dots}{further arguments passed to or from other methods.}
}

\details{The cdf is evaluated at points given by the optional argument  \code{data}. If not provided, a regular grid of length \code{ngrid} for the evaluation points is used. 

The quantiles are computed using interpolating splines on an adaptive finer grid.
}

\value{
\code{cdfMclust} returns a list of \code{x} and \code{y} values providing, respectively, the evaluation points and the estimated cdf. 

\code{quantileMclust} returns a vector of quantiles.
}

\section{References}{
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.  
}

\author{Luca Scrucca}

\seealso{
  \code{\link{densityMclust}}, 
  \code{\link{plot.densityMclust}}.
}

\examples{
x <- c(rnorm(100), rnorm(100, 3, 2))
dens <- densityMclust(x)
summary(dens, parameters = TRUE)
cdf <- cdfMclust(dens)
str(cdf)
q <- quantileMclust(dens, p = c(0.01, 0.1, 0.5, 0.9, 0.99))
cbind(quantile = q, cdf = cdfMclust(dens, q)$y)
plot(cdf, type = "l", xlab = "x", ylab = "CDF")
points(q, cdfMclust(dens, q)$y, pch = 20, col = "red3")

par(mfrow = c(2,2))
dens.waiting <- densityMclust(faithful$waiting)
plot(dens.waiting)
plot(cdfMclust(dens.waiting), type = "l", 
     xlab = dens.waiting$varname, ylab = "CDF")
dens.eruptions <- densityMclust(faithful$eruptions)
plot(dens.eruptions)
plot(cdfMclust(dens.eruptions), type = "l", 
     xlab = dens.eruptions$varname, ylab = "CDF")
par(mfrow = c(1,1))
}

\keyword{cluster}
\keyword{dplot}
