\name{emE}
\alias{emE}
\alias{emV}
\alias{emEII}
\alias{emVII}
\alias{emEEI}
\alias{emVEI}
\alias{emEVI}
\alias{emVVI}
\alias{emEEE}
\alias{emEEV}
\alias{emVEV}
\alias{emVVV}
\title{
  EM algorithm starting with E-step for a parameterized MVN mixture model.
}
\description{
  Implements the EM algorithm for a parameterized MVN mixture model,
  starting with the expectation step.
}
\usage{
emE(data, mu, sigmasq, pro, eps, tol, itmax, equalPro, warnSingular,
    Vinv, \dots)
emV(data, mu, sigmasq, pro, eps, tol, itmax, equalPro, warnSingular,
    Vinv, \dots)
emEII(data, mu, sigmasq, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emVII(data, mu, sigmasq, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emEEI(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emVEI(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emEVI(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emVVI(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emEEE(data, mu, Sigma, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emEEV(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emVEV(data, mu, decomp, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
emVVV(data, mu, sigma, pro, eps, tol, itmax, equalPro, warnSingular,
      Vinv, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables.
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the means of the
    components. 
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
  \item{Sigma}{ 
    for the equal variance model "EEE". A \emph{d} by \emph{d} matrix
    giving the common covariance for all components of the mixture
    model.
  }
  \item{sigma}{ 
    for the unconstrained variance model "VVV".
    A \emph{d} by \emph{d} by \emph{G} matrix array whose 
    \code{[,,k]}th entry is the covariance matrix for the
    \emph{k}th component of the mixture model.
  }
  \item{\dots}{
    An argument giving the variance that takes one of the following forms:
    \describe{
      \item{decomp}{for models "VVV", "EII" and "VII"; see
	\code{\link{cdens}}.} 
      \item{cholSigma}{see Sigma, for "EEE".}
      \item{cholsigma}{see sigma, for "VVV".}
      \item{sigma}{see sigma, for "VVV".}
      \item{Sigma}{see Sigma, for "EEE".}

      The form of the variance specification is the same as for the output
      for the \code{em}, \code{me}, or \code{mstep} methods for the
      specified mixture model.

      Also used to catch unused arguments from a \code{do.call} call.
    }
  }
  \item{pro}{
    Mixing proportions for the components of the mixture. 
    There should one more mixing proportion than the number of MVN components
    if the mixture model includes a  Poisson noise term.
  }
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due to
    computational singularity in covariances.  
    Smaller values of \code{eps} allow computations to
    proceed nearer to singularity. 
    The default is \code{.Mclust\$eps}.
  }
  \item{tol}{
    A scalar tolerance for relative convergence of the loglikelihood values. 
    The default is \code{.Mclust\$tol}.
  }
  \item{itmax}{
    An integer limit on the number of EM iterations. 
    The default is \code{.Mclust\$itmax}.
  }
  \item{equalPro}{
    A logical value indicating whether or not the components in the model are 
    present in equal proportions.
    The default is \code{.Mclust\$equalPro}.
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is \code{.Mclust\$warnSingular}.
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region.
    The default is determined by applying function 
    \code{hypvol} to the data.
    Used only when \code{pro} includes an additional
    mixing proportion for a noise component.
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
  \item{mu}{
    A matrix whose kth column is the
    mean of the \emph{k}th component of the mixture model.
  }
  \item{sigma}{
    An array whose \code{[,,k]}th entry is the
    the covariance matrix of the \emph{k}th component of the mixture.  
  }
  \item{pro}{
    A vector whose \emph{k}th component
    is the mixing proportion for the \emph{k}th component of the mixture model.
  }
  \item{modelName}{
    Character string identifying the model.
  }
  \item{Attributes:}{
    \itemize{
      \item \code{"info"}: Information on the iteration.
      \item \code{"warn"}: An appropriate warning if problems are
      encountered in the computations.
    }
  }
}

\details{
  This function can be used with an indirect or list call
  using \code{do.call}, allowing the output
  of e.g. \code{mstep} to be passed 
  without the need to specify individual parameters as arguments.
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
  \code{\link{em}},
  \code{\link{mstep}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}} 
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

msEst <- mstepEEE(data = irisMatrix, z = unmap(irisClass))
names(msEst)

emEEE(data = irisMatrix, mu = msEst$mu, pro = msEst$pro, cholSigma = msEst$cholSigma)
do.call("emEEE", c(list(data=irisMatrix), msEst)) ## alternative call
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
