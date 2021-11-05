\name{mclustLoglik}
\alias{mclustLoglik}
\alias{print.mclustLoglik}

\title{Log-likelihood from a table of BIC values for parameterized Gaussian mixture models}

\description{
Compute the maximal log-likelihood from a table of BIC values contained in a \code{'mclustBIC'} object as returned by function \code{\link{mclustBIC}}.
}

\usage{
mclustLoglik(object, \dots)
}

\arguments{
  \item{object}{An object of class \code{'mclustBIC'} containing the 
  BIC values as returned by a call to \code{\link{mclustBIC}}.
  }
  \item{\dots}{
    Catches unused arguments in an indirect or list call via \code{do.call}.
  }
}

\value{
An object of class \code{'mclustLoglik'} containing the maximal log-likelihood values for the Gaussian mixture models provided as input.
}

\seealso{
  \code{\link{mclustBIC}}.
}
\examples{
\donttest{
BIC <- mclustBIC(iris[,1:4])
mclustLoglik(BIC)
}
}
\keyword{cluster}
