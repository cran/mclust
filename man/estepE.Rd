\name{estepE}
\alias{estepE}
\alias{estepV}
\alias{estepEII}
\alias{estepVII}
\alias{estepEEI}
\alias{estepVEI}
\alias{estepEVI}
\alias{estepVVI}
\alias{estepEEE}
\alias{estepEEV}
\alias{estepVEV}
\alias{estepVVV}
\title{
  E-step in the EM algorithm for a parameterized MVN mixture model.
}
\description{
  Implements the expectation step in the EM algorithm for a 
  parameterized MVN mixture model.
}
\usage{
estepE(data, mu, sigmasq, pro, eps, warnSingular, Vinv, \dots)
estepV(data, mu, sigmasq, pro, eps, warnSingular, Vinv, \dots)
estepEII(data, mu, sigmasq, pro, eps, warnSingular, Vinv, \dots)
estepVII(data, mu, sigmasq, pro, eps, warnSingular, Vinv, \dots)
estepEEI(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepVEI(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepEVI(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepVVI(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepEEE(data, mu, Sigma, pro, eps, warnSingular, Vinv, \dots)
estepEEV(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepVEV(data, mu, decomp, pro, eps, warnSingular, Vinv, \dots)
estepVVV(data, mu, sigma, pro, eps, warnSingular, Vinv, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the
    means of the components.
  }
  \item{sigmasq}{
    for the one-dimensional models ("E", "V") and spherical
    models ("EII", "VII"). This is either a vector whose \emph{k}th
    component is the variance for the \emph{k}th component in the
    mixture model ("V" and "VII"), or a scalar giving the common
    variance for all components in the mixture model ("E" and "EII").
  }
  \item{decomp}{
    for the diagonal models ("EEI", "VEI", "EVI", "VVI") and some
    ellipsoidal models ("EEV", "VEV"). This is a list described in more
    detail in \code{\link{cdens}}.
  }
  \item{sigma}{
    for the unconstrained variance model "VVV" or the equal variance
    model "EEE". A \emph{d} by \emph{d} by \emph{G} matrix array whose
    \code{[,,k]}th entry is the covariance matrix for the \emph{k}th
    component of the mixture model.
  }
  \item{Sigma}{
    for the equal variance model "EEE". A \emph{d} by \emph{d} matrix
    giving the common covariance for all components of the mixture
    model.
  }
  \item{pro}{
    Mixing proportions for the components of the mixture. There should
    one more mixing proportion than the number of MVN components if the
    mixture model includes a  Poisson noise term. 
  }
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due
    to computational singularity in covariances. Smaller values of
    \code{eps} allow computations to proceed nearer to singularity. The
    default is \code{.Mclust\$eps}. 
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered. The default is
    \code{.Mclust\$warnSingular}. 
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region. The
    default is determined by applying function \code{hypvol} to the
    data. Used only when \code{pro} includes an additional mixing
    proportion for a noise component. 
  }
  \item{\dots}{
    Other arguments to describe the variance, in particular
    \code{decomp}, \code{sigma} or \code{cholsigma} for model "VVV",
    \code{decomp} for models "VII" and "EII", and \code{Sigma} or
    \code{cholSigma} for model "EEE".  \code{Sigma} is an \emph{d} by
    \emph{d} matrix giving the common covariance for all components of
    the mixture model.

    Also used to catch unused arguments from a \code{do.call} call.
  }

  }
}
\value{
  A list including the following components: 
  
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{modelName}{
    Character string identifying the model.
  }
  \item{Attribute}{
    \itemize{
      \item \code{"warn"}: An appropriate warning if problems are
      encountered in the computations.
    }
  }
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep} to be passed
  without the need to specify individual parameters as arguments. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{estep}},
  \code{\link{em}},
  \code{\link{mstep}},
  \code{\link{do.call}},
  \code{\link{mclustOptions}} 
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

msEst <- mstepEII(data = irisMatrix, z = unmap(irisClass))
names(msEst)

estepEII(data = irisMatrix, mu = msEst$mu, pro = msEst$pro, 
         sigmasq = msEst$sigmasq)
\dontrun{
do.call("estepEII", c(list(data=irisMatrix), msEst)) ## alternative call
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
