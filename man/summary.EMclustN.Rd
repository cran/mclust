\name{summary.EMclustN}
\alias{summary.EMclustN}
\title{
  summary function for EMclustN
}
\description{
  Optimal model characteristics and classification
  for \code{EMclustN} results.
}
\usage{
summary.EMclustN(object, data, G, modelNames, \dots)
}
\arguments{
  \item{object}{
    An \code{"EMclustN"} object, whch is the result of a pplying
    \code{EMclustN} to \code{data} with an initail noise estimate.
  }
  \item{data}{
    The matrix or vector of observations used to generate `object'.
  }
  \item{G}{
    A vector of integers giving the numbers of mixture components (clusters)
    over which the
    summary is to take place (\code{as.character(G)}
    must be a subset of the column names of `object').
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
  \code{\link{EMclustN}}
}

\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])

b <- apply( irisMatrix, 2, range)
n <- 450
set.seed(0)
poissonNoise <- apply(b, 2, function(x, n=n) 
                      runif(n, min = x[1]-0.1, max = x[2]+.1), n = n)
set.seed(0)
noiseInit <- sample(c(TRUE,FALSE),size=150+450,replace=TRUE,prob=c(3,1))
irisNoise <- rbind(irisMatrix, poissonNoise)

Bic <- EMclustN(data=irisNoise, noise = noiseInit)
summary(Bic, irisNoise)
summary(Bic, irisNoise, G = 0:6, modelName = c("VII", "VVI", "VVV"))
}
\keyword{cluster}
% docclass is function
