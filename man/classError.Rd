\name{classError}
\alias{classError}
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
  A list with the following two components:
  \item{misclassified}{
    The indexes of the misclassified data points in a minimum error
    mapping between the given classification and the given truth.
  }
  \item{errorRate}{
    The errorRate corresponding to a minimum error mapping 
    mapping between the given classification and the given truth.
  }
}
\details{
  If  more than one mapping between classification and truth 
  corresponds to the minimum number of classification errors,
  only one possible set of misclassified observations is returned.
}
\section{References}{
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{mapClass}},
  \code{\link{table}}
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
classError(a, b)

a <- sample(1:3, 9, replace = TRUE)
a
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
classError(a, b)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
