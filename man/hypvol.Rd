\name{hypvol}
\alias{hypvol}
\alias{pcvol}
\title{Estimation of hypervolume}
\description{
  Returns the minimum hypervolume. \code{pcvol} calculates this volume
  using principal components. \code{hypvol} takes the smaller of either
  this value, or a value calculated from the range of the data.
}

\usage{
hypvol(data, reciprocal=F)
pcvol(data, reciprocal=F)
}

\arguments{
  \item{data}{Data matrix.}
  \item{reciprocal}{A logical indicating whether the reciprocal volume
    is to be returned.}
}

\examples{
data(iris)
hypvol(iris[,1:4])
pcvol(iris[,1:4])
}

\keyword{cluster}
