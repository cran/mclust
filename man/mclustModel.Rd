\name{mclustModel}
\alias{mclustModel}
\title{
  Best model based on BIC.
}
\description{
  Determines the best model from clustering via \code{mclustBIC}
  for a given set of model parameterizations and numbers of components.
}
\usage{
mclustModel(data, BICvalues, G, modelNames, \dots)
}
\arguments{
  \item{data}{
    The matrix or vector of observations used to generate `object'.
  }
  \item{BICvalues}{
    An \code{"mclustBIC"} object, 
    which is the result of applying \code{mclustBIC} 
    to \code{data}.
  }
  \item{G}{
    A vector of integers giving the numbers of mixture components (clusters)
    from which the best model according to BIC will be selected 
   (\code{as.character(G)} must be a subset of the row names of 
     \code{BICvalues}).
    The default is to select the best model for all numbers 
    of mixture components used to obtain \code{BICvalues}.
  }
  \item{modelNames}{
    A vector of integers giving the model parameterizations
    from which the best model according to BIC will be selected 
   (\code{as.character(model)} must be a subset of the column names of 
     \code{BICvalues}).
    The default is to select the best model for parameterizations
    used to obtain \code{BICvalues}.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
 A list giving the optimal (according to BIC) parameters,
  conditional probabilities \code{z}, and loglikelihood,
  together with the associated classification and its uncertainty.

  The details of the output components are as follows:

  \item{modelName}{
    A character string denoting the model corresponding to the optimal BIC.
  }
  \item{n}{
    The number of observations in the data.
  }
  \item{d}{
    The dimension of the data.
  }
  \item{G}{
    The number of mixture components in the model corresponding to the optimal
    BIC.
  }
  \item{bic}{
    The optimal BIC value.
  }
  \item{loglik}{
    The loglikelihood corresponding to the optimal BIC.
  }
  \item{z}{
    A matrix whose \emph{[i,k]}th entry is the probability that observation
    \emph{i} in the test data belongs to the \emph{k}th class.
  }
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
  \code{\link{mclustBIC}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
mclustModel(iris[,-5], irisBIC)
mclustModel(iris[,-5], irisBIC, G = 1:6, modelNames = c("VII", "VVI", "VVV"))
}
\keyword{cluster}
% docclass is function
