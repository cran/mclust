\name{estep}
\alias{estep}
\title{
  E-step for parameterized MVN mixture models.
}
\description{
  Implements the expectation step of EM algorithm for parameterized MVN
  mixture models. 
}
\usage{
estep(modelName, data, mu, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model: \cr\cr 
    "E": equal variance  (one-dimensional) \cr
    "V": variable variance (one-dimensional) \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape \cr
    "VEI": diagonal, varying volume, equal shape\cr 
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows correspond
    to observations and columns correspond to variables. 
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the means of the
    components. 
  }
  \item{\dots}{
    Arguments for model-specific functions. Specifically:
    \itemize{
      \item An argument describing the variance (depends on the model):
      \describe{
	\item{sigmasq}{
	  for the one-dimensional models ("E", "V") and spherical models
	  ("EII", "VII"). This is either a vector whose \emph{k}th
	  component is the variance for the \emph{k}th component in the
	  mixture model ("V" and "VII"), or a scalar giving the common
	  variance for all components in the mixture model ("E" and
	  "EII").
	}
	\item{decomp}{
	  for the diagonal models ("EEI", "VEI", "EVI", "VVI") and some
	  ellipsoidal models ("EEV", "VEV"). This is a list described in
	  \code{\link{cdens}}.
	}
	\item{Sigma}{
	  for the equal variance model "EEE". A \emph{d} by \emph{d}
	  matrix giving the common covariance for all
	  components of the
	  mixture model.
	}
	\item{sigma}{
	  for the unconstrained variance model "VVV". A
	  \emph{d} by \emph{d} by \emph{G} matrix array whose
	  \code{[,,k]}th entry is the covariance matrix for
	  the \emph{k}th component of the mixture model. 
	}
	The form of the variance specification is the same as
	for the output for the \code{em}, \code{me}, or
	\code{mstep} methods for the specified mixture model. 
      }
      \item{pro}{
	Mixing proportions for the components of the
	mixture. There should one more mixing proportion than the
	number of MVN components if the mixture model includes a
	Poisson noise term. 
      }
      \item{eps}{
	A scalar tolerance for deciding when to terminate
	computations due to computational singularity in
	covariances. Smaller values of \code{eps} allow
	computations to proceed nearer to singularity. The default
	is \code{.Mclust\$eps}. 
      }
      \item{warnSingular}{
	A logical value indicating whether or not a warning
	should be issued whenever a singularity is
	encountered. The default is \code{.Mclust\$warnSingular}.  
      }
      \item{Vinv}{
	An estimate of the reciprocal hypervolume of the data region.
	The default is determined by applying function 
	\code{hypvol} to the data.
	Used only when \code{pro} includes an additional
	mixing proportion for a noise component.
      }
    }
  }
}
\value{
  A list including the following components: 
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the conditional probability
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{modelName}{
    A character string identifying the model (same as the input argument).
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
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{estepE}}, ...,
  \code{\link{estepVVV}},
  \code{\link{em}},
  \code{\link{mstep}},
  \code{\link{do.call}},
  \code{\link{mclustOptions}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

msEst <- mstep(modelName = "EII", data = irisMatrix, 
               z = unmap(irisClass))
names(msEst)

estep(modelName = msEst$modelName, data = irisMatrix,
      mu = msEst$mu, sigmasq = msEst$sigmasq, pro = msEst$pro)           
\dontrun{
do.call("estep", c(list(data = irisMatrix), msEst))   ## alternative call
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
