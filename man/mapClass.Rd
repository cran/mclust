\name{mapClass}
\alias{mapClass}

\title{Correspondence between classifications}

\description{
  Best correspondence between classes given two vectors viewed
  as alternative classifications of the same object.
}
\usage{
mapClass(a, b)
}
\arguments{
  \item{a}{
    A numeric or character vector of class labels.
  }
  \item{b}{
    A numeric or character vector of class labels.
    Must have the same length as
    \code{a}.
  }
  
}
\value{
  A list with two named elements, 
  \code{aTOb} and 
  \code{bTOa} which are themselves lists. 
  The \code{aTOb} list has a component corresponding
  to each unique element of \code{a}, which gives
  the element or elements of \code{b} 
  that result in the closest class correspondence.
  
  The \code{bTOa} list has a component corresponding
  to each unique element of \code{b}, which gives
  the element or elements of \code{a} 
  that result in the closest class correspondence.
}
\seealso{
  \code{\link{classError}},
  \code{\link{table}}
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
mapClass(a, b)
a <- sample(1:3, 9, replace = TRUE)
a
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
mapClass(a, b)
}
\keyword{cluster}
