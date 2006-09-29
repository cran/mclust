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
    The default is \code{c("E", "V")} for univariate data and
    \code{mclustOptions()\$emModelNames} for multivariate data (n > d),
    the spherical and diagonal models 
   \code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")} 
    for multivariate data (n <= d).
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
  \itemize{
  \item{hcPairs}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. For multivariate data, the default is to compute
    a hierarchical clustering tree by applying function \code{hc} with
    \code{modelName = "VVV"} to the data or a
    subset as indicated by the \code{subset} argument.
    The hierarchical clustering results are to start EM.
    For univariate data, the default is to use quantiles to start EM.
  }
  \item{subset}{
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
 A list giving the optimal (according to BIC) parameters,
  conditional probabilities \code{z}, and loglikelihood,
  together with the associated classification and its uncertainty.
  The details of the output components are as follows:
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
  \item{z}{
    A  matrix whose \emph{[i,k]}th entry is the probability that observation
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
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  
  C. Fraley and A. E. Raftery (2005).
  Bayesian regularization for normal mixture estimation and model-based
  clustering.
  Technical Report, Department of Statistics, University of Washington. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics,  University of Washington.
}
\seealso{
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustBIC}}, 
  \code{\link{mclustModelNames}},
  \code{\link{mclustOptions}}
}
\examples{
irisMclust <- Mclust(iris[,-5])
\dontrun{
 plot(irisMclust)
}
}
\keyword{cluster}
% docclass is function
