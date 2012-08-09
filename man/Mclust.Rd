\name{Mclust}
\alias{Mclust}
\alias{print.Mclust}
\title{
  Model-Based Clustering
}
\description{
  The optimal model according to BIC for EM initialized by hierarchical 
  clustering for parameterized Gaussian mixture models.
}
\usage{
Mclust(data, G=NULL, modelNames=NULL, prior=NULL, control=emControl(), 
       initialization=NULL, warn=FALSE, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the BIC is to be calculated. 
    The default is \code{G=1:9}. 
  }
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. The help file for
    \code{mclustModelNames} describes the available models.
    The default is:
    \describe{
        \item{\code{c("E", "V")}}{for univariate data}
        \item{\code{mclust.options("emModelNames")}}{for multivariate data (n > d)}
        \item{\code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}}{the spherical and diagonal models for multivariate data (n <= d)}
     }
   }
  \item{prior}{
    The default assumes no prior, but this argument allows specification of a 
    conjugate prior on the means and variances through the function 
    \code{priorControl}.
    }
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
 \item{initialization}{
    A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. For multivariate data, the default is to compute
    a hierarchical clustering tree by applying function \code{hc} with
    \code{modelName = "VVV"} to the data or a subset as indicated by the
    \code{subset} argument.
    The hierarchical clustering results are to start EM.
    For univariate data, the default is to use quantiles to start EM.
  }
  \item{\code{subset}}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.
  }
  }
  }
 \item{warn}{
    A logical value indicating whether or not certain warnings 
    (usually related to singularity) should be issued.
    The default is to suppress these warnings.
  }
  \item{\dots }{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
 An object of class \code{"Mclust"} providing the optimal (according to BIC)
 mixture model estimation.
 
 The details of the output components are as follows:
  
  \item{call}{The matched call} 

  \item{modelName}{
    A character string denoting the model at which the optimal BIC occurs.
  }
  \item{n}{
    The number of observations in the data.
  }
  \item{d}{
    The dimension of the data.
  }
  \item{G}{
    The optimal number of mixture components.
  }
  \item{BIC}{
    All BIC values.
  }
  \item{bic}{
    Optimal BIC value.
  }
  \item{loglik}{
    The loglikelihood corresponding to the optimal BIC.
  }
  \item{df}{
    The number of estimated parameters.
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
              specification. See the help file for \code{mclustVariance}
              for details.
        }
      }
  }
  \item{z}{
    A matrix whose \emph{[i,k]}th entry is the probability that observation
    \emph{i} in the test data belongs to the \emph{k}th class.
  }
  \item{classification}{
    \code{map(z)}: The classification corresponding to \code{z}.
  }
  \item{uncertainty}{
     The uncertainty associated with the classification.
  }
  \item{Attributes:}{
     The input parameters other than the data.
  }
}
\section{References}{
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 

  C. Fraley and A. E. Raftery (2005, revised 2009).
  Bayesian regularization for normal mixture estimation and model-based
  clustering.
  Technical Report, Department of Statistics, University of Washington. 

  C. Fraley and A. E. Raftery (2007).
  Bayesian regularization for normal mixture estimation and model-based
  clustering. \emph{Journal of Classification 24:155-181}.
}
\seealso{
  \code{\link{summary.Mclust}}, 
  \code{\link{plot.Mclust}}, 
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustBIC}}, 
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
mod1 = Mclust(iris[,1:4])
summary(mod1)

mod2 = Mclust(iris[,1:4], G = 1)
summary(mod2)

mod3 = Mclust(iris[,1:4], prior = priorControl())
summary(mod3)

mod4 = Mclust(iris[,1:4], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(mod4)
}
\keyword{cluster}
% docclass is function
