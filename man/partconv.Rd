\name{partconv}
\alias{partconv}

\title{Convert partitioning into numerical vector.}
\description{\code{partconv} converts a partitioning into a numerical
  vector. The second argument is used to force consecutive numbers
  (default) or not.
}
\usage{
partconv(x, consec=TRUE)
}

\arguments{
  \item{x}{Partitioning. Maybe numerical or not.}
  \item{consec}{Logical flag, whether or not to use consecutive class
    numbers.}
}
\value{
  Vector of class numbers.
}

\examples{
data(iris)
partconv(iris[,5])

cl _ sample(1:10, 25, replace=TRUE)
partconv(cl, consec=FALSE)
partconv(cl, consec=TRUE)
}

\keyword{cluster}
