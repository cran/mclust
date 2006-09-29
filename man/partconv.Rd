\name{partconv}
\alias{partconv}

\title{Numeric Encoding of a Partitioning}

\description{
  Converts a vector interpreted as a classification or partitioning 
  into a numeric vector. 
}
\usage{
partconv(x, consec=TRUE)
}

\arguments{
  \item{x}{
    A vector interpreted as a classification or partitioning. 
  }
  \item{consec}{
     Logical value indicating whether or not to consecutive class
    numbers should be used .
  }
}
\value{
  Numeric encoding of \code{x}. 
  When \code{consec = TRUE}, the distinct values in \code{x} are numbered by
  the order in which they appear.
  When \code{consec = FALSE}, each distinct value in \code{x} is numbered by
  the index corresponding to its first appearance in \code{x}.
}
\seealso{
   \code{\link{partuniq}}
}
\examples{
partconv(iris[,5])

set.seed(0)
cl <- sample(LETTERS[1:9], 25, replace=TRUE)
partconv(cl, consec=FALSE)
partconv(cl, consec=TRUE)
}

\keyword{cluster}
