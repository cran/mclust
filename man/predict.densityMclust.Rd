\name{predict.densityMclust}
\alias{predict.densityMclust}

\title{Density estimate of multivariate observations by Gaussian finite mixture modeling}

\description{Compute density estimation for multivariate observations based on Gaussian finite mixture models estimated by \code{\link{densityMclust}}.}

\usage{
  \method{predict}{densityMclust}(object, newdata, what = c("dens", "cdens"), \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'densityMclust'} resulting from a call to \code{\link{densityMclust}}.}

  \item{newdata}{a vector, a data frame or matrix giving the data. If missing the density is computed for the input data obtained from the call to \code{\link{densityMclust}}.}

  \item{what}{a character string specifying what to retrieve: \code{"dens"} returns a vector of values for the mixture density, \code{cdens} returns a matrix of component densities for each mixture component (along the columns).}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
Returns a vector or a matrix of densities evaluated at \code{newdata} depending on the argument \code{what} (see above).
}

\references{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{Mclust}}.}

\examples{
x = faithful$waiting
dens = densityMclust(x)
x0 = seq(50, 100, by = 10)
d0 = predict(dens, x0)
plot(dens)
points(x0, d0, pch = 20)
}

\keyword{multivariate}
