\name{traceW}
\alias{traceW}
\title{Compute traceW}
\description{
  Computes the trace of the sample cross product matrix.
}

\usage{traceW(x)}

\arguments{\item{x}{Data matrix}}

\examples{
data(iris)
traceW(as.matrix(iris[,1:4]))
}
\keyword{clustering}
