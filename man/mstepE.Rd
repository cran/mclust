\name{mstepE}
\alias{mstepE}
\alias{mstepV}
\alias{mstepEII}
\alias{mstepVII}
\alias{mstepEEI}
\alias{mstepVEI}
\alias{mstepEVI}
\alias{mstepVVI}
\alias{mstepEEE}
\alias{mstepEEV}
\alias{mstepVEV}
\alias{mstepVVV}
\title{
  M-step in the EM algorithm for a parameterized MVN mixture model.
}
\description{
  Maximization step in the EM algorithm for a parameterized MVN mixture model.
}
\usage{
mstepE(data, z, equalPro, noise = FALSE, \dots)
mstepV(data, z, equalPro, noise = FALSE, \dots)
mstepEII(data, z, equalPro, noise = FALSE, \dots)
mstepVII(data, z, equalPro, noise = FALSE, \dots)
mstepEEI(data, z, equalPro, noise = FALSE, eps, warnSingular, \dots)
mstepVEI(data, z, equalPro, noise = FALSE, eps, tol, itmax, warnSingular, \dots)
mstepEVI(data, z, equalPro, noise = FALSE, eps, warnSingular, \dots)
mstepVVI(data, z, equalPro, noise = FALSE, eps, warnSingular, \dots)
mstepEEE(data, z, equalPro, noise = FALSE, \dots)
mstepEEV(data, z, equalPro, noise = FALSE, eps, warnSingular, \dots)
mstepVVV(data, z, equalPro, noise = FALSE, \dots)
}
\arguments{
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
  \item{equalPro}{
    A logical value indicating whether or not the components in the model are 
    present in equal proportions.
    The default is \code{.Mclust\$equalPro}. 
  }
  \item{noise}{
    A logical value indicating whether or not the model includes a Poisson
    noise component.
    The default assumes there is no noise component.
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region.
    The default is determined by applying function 
    \code{hypvol} to the data.
    Used only when \code{noise = TRUE}.
  }
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due to
    computational singularity in covariances.
    Smaller values of \code{eps} allows computations to
    proceed nearer to singularity. 
    The default is \code{.Mclust\$eps}. 
    
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
  \item{\dots }{
    Provided to allow lists with elements other than the arguments can
    be passed in indirect or list calls with \code{do.call}.
  }
}
\value{
  A list including the following components: 
  
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
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{Attributes:}{
    \item \code{"info"} Information on the iteration.
    \item \code{"warn"} An appropriate warning if problems are
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
  \code{\link{mstep}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclustOptions}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

mstepVII(data = irisMatrix, z = unmap(irisClass))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
