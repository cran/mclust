\name{cdfMclust}
\alias{cdfMclust}

\title{
  Cumulative density function from \code{mclustDensity} estimation
}

\description{
Computes the estimated CDF from a one-dimensional density estimation at points given by the optional argument data. If not provided, a regular grid of evalutaion points is used. 
}

\usage{
cdfMclust(object, data, ngrid = 100, \dots)
}

\arguments{
  \item{object}{a \code{densityMclust} model object.}
  \item{data}{a numeric vector of evaluation points.}
  \item{ngrid}{the number of points in a regular grid to be used as evaluation points if no \code{data} are provided.}
 \item{\dots}{further arguments passed to or from other methods.}
}

\value{
Return a list of \code{x} and \code{y} values providing the evaluation points and the calculated CDF.}

%\details{}

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
plot(cdf, type = "l", xlab = "x", ylab = "CDF")

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
