\name{plot.densityMclust}
\alias{plot.densityMclust}
\title{
  Plot Univariate Mclust Density
}
\description{
  Plot method for univariate data
  associated with the \code{densityMclust} function.
}
\usage{
\method{plot}{densityMclust}(x, data = NULL, \ldots)
}
\arguments{
  \item{x}{
    Output from \code{densityMclust}.
  }
  \item{data}{
    Optional data points to be highlighted on the density curve.
  }
  \item{\dots}{
    Graphics parameters for the \code{plot} function.
  }
}
\value{
  A density curve using the mixture model and range information from
  input \code{x}. If \code{data} is supplied, the positions of the
  corresponding data points will be shown on the horizontal axis.
}
\details{
 For more flexibility in plotting, use \code{mclust1Dplot}, 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006, revised 2010).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{densityMclust}},
  \code{\link{mclust1Dplot}}
}
\examples{
\dontrun{
plot(densityMclust(rivers))
plot(densityMclust(rivers), data = rivers, xlim = c(0,1500))
}
}
\keyword{cluster}
% docclass is function
