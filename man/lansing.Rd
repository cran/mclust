\name{lansing}
\alias{lansing}
\non_function{}
\title{Maple trees in Lansing Woods}
\usage{data(lansing)}
\description{
The \code{lansing} data frame has 1217 rows and 3 columns.
The first two columns give the location, the third column the tree type.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector}
    \item{y}{a numeric vector}
    \item{species}{a factor with levels \code{hickory} and \code{maple} }
  }
}
\source{
 D.J. Gerrard, Research Bulletin No. 20, Agricultural Experimental
 Station, Michigan State University, 1969.
}
\seealso{
  \code{\link{grid1}},
  \code{\link{dens}}
}
\examples{
data(lansing)
plot(lansing[,1:2], pch=as.integer(lansing[,3]),
     col=as.integer(lansing[,3]), main="Lansing Woods tree types")
}
\keyword{datasets}
