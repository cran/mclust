\name{cdens}
\alias{cdens}
\title{
  Component Density for Parameterized MVN Mixture Models
}
\description{
  Computes component densities for observations in parameterized MVN mixture
  models.
}
\usage{
cdens(modelName, data, mu, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. Possible models: \cr\cr
    "E" for spherical, equal variance (one-dimensional) \cr
    "V" for spherical, variable variance (one-dimensional) \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume, equal shape \cr
    "VEI": diagonal, varying volume, equal shape \cr
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume, varying shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation \cr\cr
    For fitting a single Gaussian:\cr\cr    
    "X": one-dimensional \cr
    "XII": spherical \cr
    "XXI": diagonal \cr
    "XXX": ellipsoidal
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the means of the components.
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
	  ellipsoidal models ("EEV", "VEV"). This is a list with the following
	  components:
	  \describe{
	    \item{d}{
	      The dimension of the data.
	    }
	    \item{G}{
	      The number of components in the mixture model.
	    }
	    \item{scale}{
	      Either a \emph{G}-vector giving the scale of the covariance (the
	      \emph{d}th root of its determinant) for each component in the
	      mixture model, or a single numeric value if the scale is the
	      same for each component.
	    }
	    \item{shape}{
	      Either a \emph{G} by \emph{d} matrix in which the \emph{k}th
	      column is the shape of the covariance matrix (normalized to have
	      determinant 1) for the \emph{k}th component, or a
	      \emph{d}-vector giving a common shape for all components.
	    }
	    \item{orientation}{
	      Either a \emph{d} by \emph{d} by \emph{G} array whose
	      \code{[,,k]}th entry is the orthonomal matrix of eigenvectors of
	      the covariance matrix of the \emph{k}th component, or a \emph{d}
	      by \emph{d} orthonormal matrix if the mixture components have a
	      common orientation. The \code{orientation} component of
	      \code{decomp} can be omitted in spherical and diagonal models,
	      for which the principal components are parallel to the
	      coordinate axes so that the orientation matrix is the identity.
	    }
	  }
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
      \item \code{eps}:
      A scalar tolerance for deciding when to terminate computations
      due to computational singularity in covariances. Smaller values
      of \code{eps} allow computations to proceed nearer to
      singularity. The default is \code{.Mclust\$eps}. 
      
      For those models with iterative M-step ("VEI", "VEV"), two
      values can be entered for \code{eps}, in which case the second
      value is used for determining singularity in the M-step. 
      
      \item \code{warnSingular}:
      A logical value indicating whether or not a warning should be issued
      whenever a singularity is encountered.
      The default is \code{.Mclust\$warnSingular}.
    }      
  }
}
\value{
  A numeric matrix whose \code{[i,j]}th entry is the 
  density of observation \emph{i} in component \emph{j}. 
  The densities are not scaled by mixing proportions.
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
  \code{\link{cdensE}}, \ldots,
  \code{\link{cdensVVV}},
  \code{\link{dens}},
  \code{\link{EMclust}},
  \code{\link{mstep}},
  \code{\link{mclustDAtrain}},
  \code{\link{mclustDAtest}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{
n <- 100 ## create artificial data

set.seed(0)
x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))
clPairs(x, cl = xclass, sym = c("1","2")) ## display the data

set.seed(0)
I <- sample(1:(2*n)) ## random ordering of the data
x <- x[I, ]
xclass <- xclass[I]

odd <- seq(1, 2*n, by = 2)
oddBic <- EMclust(x[odd, ]) 
oddSumry <- summary(oddBic, x[odd, ]) ## best parameter estimates
names(oddSumry)

even <- odd + 1
temp <- cdens(modelName = oddSumry$modelName, data = x[even, ],
              mu = oddSumry$mu, decomp = oddSumry$decomp)
cbind(class = xclass[even], temp)

## alternative call 

temp <- do.call( "cdens", c(list(data = x[even, ]), oddSumry))
cbind(class = xclass[even], temp)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
