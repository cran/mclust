\name{em}
\alias{em}
\title{
  EM algorithm starting with E-step for parameterized MVN mixture models.
}
\description{
  Implements the EM algorithm for parameterized MVN mixture models,
  starting with the expectation step.
}
\usage{
em(modelName, data, mu, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model: \cr\cr
    "E": equal variance  (one-dimensional) \cr
    "V": variable variance (one-dimensional) \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape\cr 
    "VEI": diagonal, varying volume, equal shape\cr 
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation\cr 
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation
  }
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
  \item{\dots}{
    Arguments for model-specific em functions. Specifically:
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
	  ellipsoidal models ("EEV", "VEV"). For a description,
	  see \code{\link{cdens}}.
	}
	\item{Sigma}{
	  for the equal variance model "EEE". A \emph{d} by \emph{d}
	  matrix giving the common covariance for all components of the
	  mixture model.
	}
	\item{sigma}{
	  for the unconstrained variance model "VVV". A \emph{d} by
	  \emph{d} by \emph{G} matrix array whose \code{[,,k]}th entry
	  is the covariance matrix for the \emph{k}th component of the
	  mixture model.
	}
	The form of the variance specification is the same as for the
	output for the \code{em}, \code{me}, or \code{mstep} methods for
	the specified mixture model. 
      }
      
      \item \code{pro}:
      Mixing proportions for the components of the mixture. There
      should one more mixing proportion than the number of MVN
      components if the mixture model includes a  Poisson noise term. 
      \item \code{eps}:
      A scalar tolerance for deciding when to terminate computations
      due to computational singularity in covariances. Smaller values
      of \code{eps} allow computations to proceed nearer to
      singularity. The default is \code{.Mclust\$eps}. 
      
      For those models with iterative M-step ("VEI", "VEV"), two
      values can be entered for \code{eps}, in which case the second
      value is used for determining singularity in the M-step. 

      \item \code{tol}:
      A scalar tolerance for relative convergence of the loglikelihood. 
      The default is \code{.Mclust\$tol}.
      
      For those models with iterative M-step ("VEI", "VEV"), two
      values can be entered for \code{tol}, in which case the second
      value governs parameter convergence in the M-step. 

      \item \code{itmax}:
      An integer limit on the number of EM iterations. 
      The default is \code{.Mclust\$itmax}.
      
      For those models with iterative M-step ("VEI", "VEV"), two
      values can be entered for \code{itmax}, in which case the second
      value is an upper limit on the number of iterations in the
      M-step. 

      \item \code{equalPro}:
      Logical variable indicating whether or not the mixing proportions are
      equal in the model. The default is \code{.Mclust\$equalPro}.
      
      \item \code{warnSingular}:
      A logical value indicating whether or not a warning should be issued
      whenever a singularity is encountered.
      The default is \code{.Mclust\$warnSingular}.
      
      \item \code{Vinv}:
      An estimate of the reciprocal hypervolume of the data region.
      The default is determined by applying function 
      \code{hypvol} to the data.
      Used only when \code{pro} includes an additional
      mixing proportion for a noise component.
    }
  }
}

\value{
  A list including the following components: 
  
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the conditional probability of
    the \emph{i}th observation belonging to the \emph{k}th component of the
    mixture.   
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{mu}{
    A matrix whose kth column is the
    mean of the \emph{k}th component of the mixture model.
  }
  \item{sigma}{
    For multidimensional models, a three dimensional array 
    in which the \code{[,,k]}th entry gives the
    the covariance for the \emph{k}th group in the best model. <br>
    For one-dimensional models, either a scalar giving a common variance for
    the groups or a vector whose entries are the variances for each group
    in the best model.
  }
  \item{pro}{
    A vector whose \emph{k}th component
    is the mixing proportion for the \emph{k}th component of the mixture model.
  }
  \item{modelName}{
    A character string identifying the model (same as the input argument).
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
  \code{\link{emE}}, \ldots,
  \code{\link{emVVV}},
  \code{\link{estep}},
  \code{\link{me}},
  \code{\link{mstep}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]
 
msEst <- mstep(modelName = "EEE", data = irisMatrix, 
               z = unmap(irisClass))
names(msEst)

em(modelName = msEst$modelName, data = irisMatrix,
   mu = msEst$mu, Sigma = msEst$Sigma, pro = msEst$pro)
do.call("em", c(list(data = irisMatrix), msEst))   ## alternative call
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
