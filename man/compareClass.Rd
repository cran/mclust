\name{compareClass}
\alias{compareClass}
\title{
  Compare classifications.
}
\description{
  Compare classifications via the normalized 
  variation of information criterion.
}
\usage{
compareClass(a, b)
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
  The variation of information criterion (Meila 2002) for \code{a} and
  \code{b} divided by the log of the length of the sequences so that it
  falls in \emph{[0,1]}.
}
\section{References}{
  Marina Meila (2002).
  Comparing clusterings.
  Technical Report no. 418, Department of Statistics, University
  of Washington.
  
  See \url{http://www.stat.washington.edu/www/research/reports}.
}
\seealso{
  \code{\link{mapClass}},
  \code{\link{classError}},
  \code{\link{table}}
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
compareClass(a, b)
a <- sample(1:3, 9, replace = TRUE)
a
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
compareClass(a, b)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
