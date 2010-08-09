\name{summary.mclustBIC}
\alias{summary.mclustBIC}
\alias{print.summary.mclustBIC}
\alias{summary.mclustBIC}
\alias{summaryMclustBIC}
\alias{summaryMclustBICn}
\alias{printSummaryMclustBIC}
\alias{printSummaryMclustBICn}
\title{
  Summary Function for model-based clustering.
}
\description{
  Optimal model characteristics and classification for model-based
  clustering via \code{mclustBIC}.
}
\usage{
\method{summary}{mclustBIC}(object, data, G, modelNames, \dots)
}
\arguments{
  \item{object}{
    An \code{"mclustBIC"} object, 
    which is the result of applying \code{mclustBIC} 
    to \code{data}.
  }
  \item{data}{
    The matrix or vector of observations used to generate `object'.
  }
 \item{G}{
    A vector of integers giving the numbers of mixture components (clusters)
    from which the best model according to BIC will be selected
   (\code{as.character(G)} must be a subset of the row names of
     \code{object}).
    The default is to select the best model for all numbers
    of mixture components used to obtain \code{object}.
  }
  \item{modelNames}{
    A vector of integers giving the model parameterizations
    from which the best model according to BIC will be selected
   (\code{as.character(model)} must be a subset of the column names of
     \code{object}).
    The default is to select the best model for parameterizations
    used to obtain \code{object}.
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
  \item{parameters}{
     A list with the following components:
     \describe{
        \item{pro}{
              A vector whose \emph{k}th component is the mixing proportion for
              the \emph{k}th component of the mixture model.
              If missing, equal proportions are assumed.
        }
        \item{mean}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{variance}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{mclustVariance}
              for details.
        }
      }
   }   
  \item{z}{
    A matrix whose \emph{[i,k]}th entry is the probability that observation
    \emph{i} in the data belongs to the \emph{k}th class.
  }
  \item{classification}{
    \code{map(z)}: The classification corresponding to \code{z}.
  }
  \item{uncertainty}{
     The uncertainty associated with the classification.
  }
  \item{Attributes:}{
    \itemize{
      \item \code{"bestBICvalues"}
        Some of the best bic values for the analysis.
      \item \code{"prior"}
        The prior as specified in the input.
      \item \code{"control"}
        The control parameters for EM as specified in the input.
      \item \code{"initialization"}
        The parameters used to initial EM for computing the maximum
        likelihood values used to obtain the BIC.
    }
  }
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}.

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}

\seealso{
  \code{\link{mclustBIC}}
  \code{\link{mclustModel}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])
summary(irisBIC, iris[,-5], G = 1:6, modelNames = c("VII", "VVI", "VVV"))
}
\keyword{cluster}
% docclass is function
