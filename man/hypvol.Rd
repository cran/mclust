\name{hypvol}
\alias{hypvol}
\title{
  Aproximate Hypervolume for Multivariate Data
}
\description{
  Computes a simple approximation to the hypervolume of a multivariate
  data set.
}
\usage{
hypvol(data, reciprocal=FALSE)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{reciprocal}{
    A logical variable indicating whether or not the reciprocal
    hypervolume is desired rather than the hypervolume itself. The
    default is to return the approximate hypervolume.  
  }
}
\value{
  Computes the hypervolume by two methods: simple variable bounds
  and principal components, and returns the minimum value.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
hypvol(irisMatrix)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
