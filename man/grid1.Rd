\name{grid1}
\alias{grid1}
\alias{grid2}
\title{Generate grid points}
\description{
Generate grid points in one or two dimensions.
}
\usage{
grid1(n, range = c(0, 1), edge = TRUE)
grid2(x, y)
}
\arguments{
  \item{n}{Number of grid points.}
  \item{range}{Range of grid points.}
  \item{edge}{Logical: include edges or not?}
  \item{x, y}{Vectors.}
}
\value{
The value returned is simple: \code{grid1} generates a
vector; \code{grid2} generates a matrix.
}
\author{C. Fraley}
\seealso{
  \code{\link{lansing}},
  \code{\link{dens}}
}
\examples{
data(lansing)
maples <- lansing[as.character(lansing[,"species"]) == "maple", -3]
maplesBIC <- EMclust(maples)
maplesModel <- summary(maplesBIC, maples)
x <- grid1(100, range=c(0,1))
y <- x
xyDens <- do.call("dens", c(list(data=grid2(x, y)), maplesModel))
xyDens <- matrix(xyDens, ncol=100)
contour(xyDens)
points(maples, cex=.2, col="red")

image(xyDens)
points(maples, cex=.5)
}
\keyword{cluster}

