\name{compClass}
\alias{compClass}
\title{
  Compare classifications having the same number of groups.
}
\description{
  Error rate and best correspondence between two 
  classifications having the same number of groups.
}
\usage{
compClass(a, b)
}
\arguments{
  \item{a}{
    A numeric or character vector of class labels.
  }
  \item{b}{
    A numeric or character vector of class labels. Must have the same
    length and the same number of class labels as \code{a}. 
  }
  
}
\value{
  A list consisting of the following components:
  
  \item{error}{
    The minimum error rate between the two classifications, 
    taken over all possible permutations. 
  }
  \item{map}{
    A matrix giving a correspondence between the unique elements of
    \code{a} and \code{b} that gives the minimum error rate. (Note:
    there may be more than one mapping that achieves the minimum error
    rate.) 
  }
}
\section{References}{
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
compClass(a, b)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.


