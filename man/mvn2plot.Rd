\name{mvn2plot}
\alias{mvn2plot}
\title{Plot confidence ellipses}
\description{
  Function plots one standard deviation of an MVN mixture classification.
  Not really meant to be called by users...
}
\usage{
mvn2plot(mu, sigma, k = 15, alone = FALSE)
}
\arguments{
  \item{mu}{Mean.}
  \item{sigma}{Covariance matrix.}
  \item{k}{Number of segments for plotting segments of
    ellipsoids. Default: 15.}
  \item{alone}{Logical: if \code{TRUE} then the coordinates of the plot
    are determined by the size of the ellipse.}
}
\seealso{\code{\link{coordProj}}}
\keyword{cluster}
