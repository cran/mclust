\name{classError}
\alias{classError}
\alias{classErrors}
\title{
  Classification error.
}
\description{
  Error for a given classification relative to a known truth. 
  Location of errors in a given classification relative to a known truth.
}
\usage{
classError(classification, truth)
}
\arguments{
  \item{classification}{
    A numeric or character vector of class labels.
  }
  \item{truth}{
    A numeric or character vector of class labels.
    Must have the same length as
    \code{classification}.
  }
  
}
\value{
  \code{classError} gives the fraction of elements 
  misclassified for \code{classification}
  relative to \code{truth}. 
  \code{classErrors} is a logical vector of the
  same length as \code{classification} and
  \code{truth} which gives the location
  of misclassified elements in \code{classification}
  relative to \code{truth}. 
  
}
\details{
  \code{classErrors} will only return one possibility
  if more than one mapping between classification and truth results in the
  minimum error.
}
\seealso{
  \code{\link{compareClass}},
  \code{\link{mapClass}},
  \code{\link{table}}
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
classError(a, b)
classErrors(a, b)

a <- sample(1:3, 9, replace = TRUE)
a
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
classError(a, b)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
