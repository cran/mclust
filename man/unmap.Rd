\name{unmap}
\alias{unmap}
\title{
  Indicator Variables given Classification
}
\description{
  Converts a classification into a matrix of indicator variables.
}
\usage{
  unmap(classification, noise, \dots)
}
\arguments{
  \item{classification}{
    A numeric or character vector. Typically the distinct entries of this
    vector would represent a classification of observations in a data set.
  }
  \item{noise}{
    A single numeric or character value used to indicate the value of
    \code{classification} corresponding to noise.
  }
 \item{\dots }{
     Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  An \emph{n} by \emph{m} matrix of \emph{(0,1)} indicator variables,
  where \emph{n} is the length of \code{classification} and \emph{m} is
  the number of unique values or symbols in  \code{classification}. 
  Columns are labeled by the unique values in \code{classification}, 
  and the \code{[i,j]}th entry is \emph{1} if \code{classification[i]} 
  is the \emph{j}th unique value or symbol in sorted order 
  \code{classification}. 
  If a \code{noise} value of symbol is designated, the corresponding indicator 
  variables are relocated to the last column of the matrix. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{map}},
  \code{\link{estep}},
  \code{\link{me}}
}
\examples{
z <- unmap(iris[,5])
z[1:5, ]
  
emEst <- me(modelName = "VVV", data = iris[,-5], z = z)
emEst$z[1:5,]
  
map(emEst$z)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
