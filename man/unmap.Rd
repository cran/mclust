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
    A single numeric or character value used to indicate observations 
    corresponding to noise.
  }
  \item{\dots }{
    Provided to allow lists with elements other than the arguments can
    be passed in indirect or list calls with \code{do.call}.
  }
}
\value{
  An \emph{n} by \emph{m} matrix of \emph{(0,1)} indicator variables,
  where \emph{n} is the length of \code{classification} and \emph{m} is
  the number of unique values or symbols in
  \code{classification}. Columns are labeled by the unique values in
  \code{classification}, and the \code{[i,j]}th entry is \emph{1} if
  \code{classification[i]} is the \emph{j}th unique value or symbol in
  order of appearance in the \code{classification}. If a \code{noise}
  value of symbol is designated, the corresponding indicator variables
  are located in the last column of the matrix. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis.
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{map}},
  \code{\link{estep}},
  \code{\link{me}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]
  
z <- unmap(irisClass)
z
  
emEst <- me(modelName = "VVV", data = irisMatrix, z = z)
emEst$z
  
map(emEst$z)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
