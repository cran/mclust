\name{meE}
\alias{meE}
\alias{meV}
\alias{meEII}
\alias{meVII}
\alias{meEEI}
\alias{meVEI}
\alias{meEVI}
\alias{meVVI}
\alias{meEEE}
\alias{meEEV}
\alias{meVEV}
\alias{meVVV}
\title{
  EM algorithm starting with M-step for a parameterized Gaussian mixture model.
}
\description{
  Implements the EM algorithm for a parameterized Gaussian mixture model,
  starting with the maximization step.
}
\usage{
meE(data, z, prior=NULL, control=emControl(), 
    Vinv=NULL, warn=NULL, \dots)
meV(data, z, prior=NULL, control=emControl(),
    Vinv=NULL, warn=NULL, \dots)
meEII(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meVII(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meEEI(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meVEI(data, z, prior=NULL, control=emControl(),
     Vinv=NULL, warn=NULL, \dots)
meEVI(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meVVI(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meEEE(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meEEV(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meVEV(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
meVVV(data, z, prior=NULL, control=emControl(),
      Vinv=NULL, warn=NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the ith observation belonging to
    the \emph{k}th component of the mixture.  
  }
 \item{prior}{
    Specification of a conjugate prior on the means and variances.
    The default assumes no prior.                                                 }
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
 \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region, when the
    model is to include a noise term. Set to a negative value or zero if
    a noise term is desired, but an estimate is unavailable --- in that
    case function \code{hypvol} will be used to obtain the estimate.
    The default is not to assume a noise term in the model through the
    setting \code{Vinv=NULL}.
  }
  \item{warn}{
    A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when the
    estimation fails. The default is given by \code{mclust.options("warn")}.
  }
 \item{\dots }{
   Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A list including the following components: 
 \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{parameters}{
     \describe{
        \item{\code{pro}}{
              A vector whose \emph{k}th component is the mixing proportion for 
              the \emph{k}th component of the mixture model.
              If the model includes a Poisson term for noise, there 
              should be one more mixing proportion than the number 
              of Gaussian components.
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
  \item{loglik}{
    The log likelihood for the data in the mixture model. 
  }
  \item{Attributes:}{
      \code{"info"} Information on the iteration.\cr
      \code{"WARNING"} An appropriate warning if problems are encountered 
      in the computations.\cr
  }
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  
  C. Fraley and A. E. Raftery (2005).
  Bayesian regularization for normal mixture estimation and model-based
  clustering.
  Technical Report, Department of Statistics, University of Washington.

  C. Fraley and A. E. Raftery (2007).
  Bayesian regularization for normal mixture estimation and model-based
  clustering. \emph{Journal of Classification 24:155-181}.

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{em}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclust.options}}
}
\examples{
meVVV(data = iris[,-5], z = unmap(iris[,5]))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
