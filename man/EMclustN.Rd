\name{EMclustN}
\alias{EMclustN}
\alias{print.EMclustN}
\alias{plot.EMclustN}
\title{
  BIC for Model-Based Clustering with Poisson Noise
}
\description{
  BIC for EM initialized by hierarchical clustering
  for parameterized Gaussian mixture models with Poisson noise.
}
\usage{
EMclustN(data, G, emModelNames, noise, hcPairs, eps, tol, itmax, 
         equalPro, warnSingular=FALSE, Vinv, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{G}{
    An integer vector specifying the numbers of MVN (Gaussian) mixture
    components (clusters) for which the BIC is to be calculated. The
    default is \code{0:9} where \code{0} indicates only a noise
    component. 
  }
  \item{emModelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. Possible models: \cr\cr    
    "E" for spherical, equal variance (one-dimensional) \cr
    "V" for spherical, variable variance (one-dimensional) \cr
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
    The default is \code{.Mclust\$emModelNames}.
  }
  \item{noise}{
    A logical vector indicating whether or not observations are initially 
    estimated to noise in the data.
    If there is no noise \code{EMclust}
    should be use rather than \code{EMclustN}.
  }
  \item{hcPairs}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. The default is to compute a hierarchical
    clustering tree by applying function \code{hc} with
    \code{modelName = .Mclust\$hcModelName[1]} to univariate data and
    \code{modelName = .Mclust\$hcModelName[2]} to multivariate data or a
    subset as indicated by the \code{subset} argument. The hierarchical
    clustering results are used as starting values for EM.  
  }
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due to
    computational singularity in covariances. Smaller values of
    \code{eps} allow computations to proceed nearer to singularity. The
    default is \code{.Mclust\$eps}. 
  }
  \item{tol}{
    A scalar tolerance for relative convergence of the loglikelihood. 
    The default is \code{.Mclust\$tol}.
  }
  \item{itmax}{
    An integer limit on the number of EM iterations. 
    The default is \code{.Mclust\$itmax}.
  }
  \item{equalPro}{
    Logical variable indicating whether or not the mixing proportions are
    equal in the model. The default is \code{.Mclust\$equalPro}.
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region.
    The default is determined by applying function 
    \code{hypvol} to the data.
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is \code{warnSingular=FALSE}.
  }
  \item{\dots }{
    Provided to allow lists with elements other than the arguments can
    be passed in indirect or list calls with \code{do.call}.
  }
}
\value{
  Bayesian Information Criterion for the specified mixture models
  numbers of clusters. Auxiliary information returned as attributes.
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
  \code{\link{summary.EMclustN}}, 
  \code{\link{EMclust}}, 
  \code{\link{hc}},
  \code{\link{me}},
  \code{\link{mclustOptions}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

b <- apply( irisMatrix, 2, range)
n <- 450
set.seed(0)
poissonNoise <- apply(b, 2, function(x, n=n) 
                      runif(n, min = x[1]-0.1, max = x[2]+.1), n = n)
set.seed(0)
noiseInit <- sample(c(TRUE,FALSE),size=150+450,replace=TRUE,prob=c(3,1))
Bic <-  EMclustN(data=rbind(irisMatrix, poissonNoise), noise = noiseInit)
Bic
plot(Bic)
}
\keyword{cluster}
% docclass is function
