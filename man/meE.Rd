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
  EM algorithm starting with M-step for a parameterized MVN mixture model.
}
\description{
  Implements the EM algorithm for a parameterized MVN mixture model,
  starting with the maximization step.
}
\usage{
meE(data, z, eps, tol, itmax, equalPro, warnSingular, 
    noise = FALSE, Vinv)
meV(data, z, eps, tol, itmax, equalPro, warnSingular, 
    noise = FALSE, Vinv)
meEII(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meVII(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meEEI(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meVEI(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meEVI(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meVVI(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meEEE(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meEEV(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meVEV(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
meVVV(data, z, eps, tol, itmax, equalPro, warnSingular, 
      noise = FALSE, Vinv)
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
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due to
    computational singularity in covariances.  
    Smaller values of \code{eps} allows computations to
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
    Logical variable indicating whether or not the mixing proportions are
    equal in the model.
    The default is \code{.Mclust\$equalPro}.
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is \code{.Mclust\$warnSingular}.
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
}
\value{
  A list including the following components: 
  \item{mu}{
    A matrix whose kth column is the mean of the \emph{k}th
    component of the mixture model.
  }
  \item{sigma}{
    For multidimensional models, a three dimensional array 
    in which the \code{[,,k]}th entry gives the
    the covariance for the \emph{k}th group in the best model. <br>
    For one-dimensional models, either a scalar giving a common variance for
    the groups or a vector whose entries are the variances for each group
    in the best model.
  }
  \item{pro}{ A vector whose \emph{k}th component
    is the mixing proportion for the \emph{k}th component of the mixture model.
  }
  \item{z}{ 
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{modelName}{
    Character string identifying the model.
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{Attributes:}{ 
    The return value also has the following attributes:
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
  \code{\link{em}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclustOptions}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

meVVV(data = irisMatrix, z = unmap(irisClass))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
