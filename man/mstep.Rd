\name{mstep}
\alias{mstep}
\title{
  M-step in the EM algorithm for parameterized MVN mixture models.
}
\description{
  Maximization step in the EM algorithm for parameterized MVN mixture models.
}
\usage{
mstep(modelName, data, z, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model: \cr\cr
    "E": equal variance  (one-dimensional) \cr
    "V": variable variance (one-dimensional) \cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape \cr
    "VEI": diagonal, varying volume, equal shape \cr
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the ith observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{\dots}{
    Any number of the following:
    \describe{
      \item{equalPro}{
	A logical value indicating whether or not the components in the
	model are present in equal proportions. The default is
	\code{.Mclust\$equalPro}.  
      }
      \item{noise}{
	A logical value indicating whether or not the model includes a Poisson
	noise component. The default assumes there is no noise component.
      }
      \item{eps}{
	A scalar tolerance for deciding when to terminate computations due
	to computational singularity in covariances. Smaller values of
	\code{eps} allows computations to proceed nearer to singularity. The
	default is \code{.Mclust\$eps}.  
	
	Not used for models "EII", "VII", "EEE", "VVV".
      }
      \item{tol}{
	For models with iterative M-step ("VEI", "VEE", "VVE", "VEV"), 
	a scalar tolerance for relative convergence of the parameters. 
	The default is \code{.Mclust\$tol}. 
      }
      \item{itmax}{
	For models with iterative M-step ("VEI", "VEE", "VVE", "VEV"), 
	an integer limit on the number of EM iterations. 
	The default is \code{.Mclust\$itmax}. 
      }
      \item{warnSingular}{
	A logical value indicating whether or not a warning should be issued
	whenever a singularity is encountered.
	The default is \code{.Mclust\$warnSingular}.
	
	Not used for models "EII", "VII", "EEE", "VVV".
      }
    }
  }
}
\value{
  A list including the following components: 
  
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
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{Attributes:}{
    \item \code{"info"}: Information on the iteration.
    \item \code{"warn"}: An appropriate warning if problems are
    encountered in the computations. 
  }
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
  \code{\link{mstepE}}, \ldots,
  \code{\link{mstepVVV}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclustOptions}}.
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]
 
mstep(modelName = "VII", data = irisMatrix, z = unmap(irisClass))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
