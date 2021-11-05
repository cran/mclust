\name{dupPartition}
\alias{dupPartition}

\title{Partition the data by grouping together duplicated data}

\description{
Duplicated data are grouped together to form a basic partition that can be used to start hierarchical agglomeration.
}
\usage{
dupPartition(data)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    If a matrix or data frame, rows correspond to observations (\eqn{n}) and
    columns correspond to variables (\eqn{d}).
  }
}
\value{
  A vector of indices indicating the partition. 
}
\seealso{
  \code{\link{hc}}
}
\examples{
\donttest{
dupPartition(iris[,1:4])
dupPartition(iris)
dupPartition(iris$Species)
}
}
\keyword{cluster}
