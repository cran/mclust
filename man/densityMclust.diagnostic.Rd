\name{densityMclust.diagnostic}
\alias{densityMclust.diagnostic}

\title{Diagnostic plots for \code{mclustDensity} estimation}

\description{
Diagnostic plots for density estimation. Only available for the one-dimensional case.
}

\usage{
densityMclust.diagnostic(object, type = c("cdf", "qq"), 
                         col = c("black", "black"), 
                         lwd = c(2,1), lty = c(1,1), 
                         legend = TRUE, grid = TRUE, 
                         \dots)
}

\arguments{
  \item{object}{An object of class \code{'mclustDensity'} obtained from a call to \code{\link{densityMclust}} function.}

  \item{type}{The type of graph requested:
    \describe{
      \item{\code{"cdf"} =}{a plot of the estimated CDF versus the empirical distribution function.}
      \item{\code{"qq"} =}{a Q-Q plot of sample quantiles versus the quantiles obtained from the inverse of the estimated cdf.}
    }
  }
  
  \item{col}{A pair of values for the color to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{lwd}{A pair of values for the line width to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{lty}{A pair of values for the line type to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{legend}{A logical indicating if a legend must be added to the plot of fitted CDF vs the empirical CDF.}

  \item{grid}{A logical indicating if a \code{\link{grid}} should be added to the plot.}

  \item{\dots}{Additional arguments.}
}

\details{
The two diagnostic plots for density estimation in the one-dimensional case are discussed in Loader (1999, pp- 87-90).

}

% \value{}

\references{
  Loader C. (1999), Local Regression and Likelihood. New York, Springer.

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
\donttest{
x <- faithful$waiting
dens <- densityMclust(x, plot = FALSE)
plot(dens, x, what = "diagnostic")
# or
densityMclust.diagnostic(dens, type = "cdf")
densityMclust.diagnostic(dens, type = "qq")
}
}
\keyword{cluster}
\keyword{dplot}
