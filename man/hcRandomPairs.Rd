\name{hcRandomPairs}
\alias{hcRandomPairs}
\alias{randomPairs}

\title{Random hierarchical structure}

\description{Create a hierarchical structure using a random hierarchical partition of the data.}

\usage{
hcRandomPairs(data, seed = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric matrix or data frame of observations.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{seed}{
    Optional single value, interpreted as an integer, specifying the seed for random partition.
  }
  \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}

\value{
  A numeric two-column matrix in which the \emph{i}th row gives the minimum 
  index for observations in each of the two clusters merged at the
  \emph{i}th stage of a random agglomerative hierarchical clustering.
}

\seealso{
  \code{\link{hc}},
  \code{\link{hclass}}
  \code{\link{hcVVV}}
}

\examples{
data <- iris[,1:4]
randPairs <- hcRandomPairs(data)
str(randPairs)
# start model-based clustering from a random partition
mod <- Mclust(data, initialization = list(hcPairs = randPairs))
summary(mod)
}

\keyword{cluster}
