\name{mclustModel}
\alias{mclustModel}
\title{
  Best model based on BIC
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
    An \code{'mclustBIC'} object, 
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
  conditional probabilities \code{z}, and log-likelihood,
  together with the associated classification and its uncertainty.

  The details of the output components are as follows:
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
   \item{n}{
    The number of observations in the data.
  }
  \item{d}{
    The dimension of the data.
  }
  \item{G}{
    The number of components in the Gaussian mixture model corresponding
    to the optimal BIC.
  }
  \item{bic}{
    The optimal BIC value.
  }
  \item{loglik}{
    The log-likelihood corresponding to the optimal BIC.
  }
 \item{parameters}{
     A list with the following components:
     \describe{
        \item{\code{pro}}{
              A vector whose \emph{k}th component is the mixing proportion for
              the \emph{k}th component of the mixture model.
              If missing, equal proportions are assumed.
        }
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}}
              for details.
            }
	\item{\code{Vinv}}{
              The estimate of the reciprocal hypervolume of the data region
              used in the computation when the input indicates the
              addition of a noise component to the model.
        }
      }
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

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
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
