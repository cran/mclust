\name{summary.EMclust}
\alias{summary.EMclust}
\alias{print.summary.EMclust}
\title{
  Summary function for EMclust
}
\description{
  Optimal model characteristics and classification
  for \code{EMclust} results.
}
\usage{
summary.EMclust(object, data, G, modelNames, \dots)
}
\arguments{
  \item{object}{
    An \code{"EMclust"} object, 
    which is the result of applying \code{EMclust} 
    to \code{data}.
  }
  \item{data}{
    The matrix or vector of observations used to generate `object'.
  }
  \item{G}{
    A vector of integers giving the numbers of mixture components (clusters)
    over which the
    summary is to take place (\code{as.character(G)}
    must be a subset of the column names of \code{object}).
    The default is to summarize over all of the numbers of mixture components 
    used in the original analysis.
  }
  \item{modelNames}{
    A vector of character strings denoting the models over which the
    summary is to take place (must be a subset of the row names of `object').
    The default is to summarize over all models used in the original
    analysis.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
  A list giving the optimal (according to BIC) parameters, 
  conditional probabilities \code{z}, and loglikelihood,
  together with the associated classification and its uncertainty.
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
  \code{\link{EMclust}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])

irisBic <- EMclust(irisMatrix)
summary(irisBic, irisMatrix)
summary(irisBic, irisMatrix, G = 1:6, modelName = c("VII", "VVI", "VVV"))
}
\keyword{cluster}
% docclass is function
