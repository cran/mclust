\name{plot.emclust}
\alias{plot.emclust1}

\title{Plot BIC values}
\description{

}
\usage{
plot.emclust(x, ...)
plot.emclust1(x, ...)
}
\arguments{
  \item{x}{\code{emclust} object containing BIC values.}
  \item{\dots}{Further arguments to the plot function.}
}

\seealso{\code{\link{emclust}}, \code{\link{emclust1}}}

\examples{
data(iris)
bicvals _ emclust(iris[,1:4], nclus=1:3, modelid=c("VVV","EEV","VEV"))
plot(bicvals, col=c("red","cyan","purple"), lty=1)
}
\keyword{cluster}
