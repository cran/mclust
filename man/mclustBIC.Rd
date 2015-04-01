\name{mclustBIC}
\alias{mclustBIC}
\alias{EMclust}
\alias{print.mclustBIC}

\title{BIC for Model-Based Clustering}

\description{
BIC for parameterized Gaussian mixture models fitted by EM algorithm initialized by model-based hierarchical clustering.
}

\usage{
mclustBIC(data, G = NULL, modelNames = NULL, 
          prior = NULL, control = emControl(), 
          initialization = list(hcPairs = NULL, 
                                subset = NULL, 
                                noise = NULL), 
          Vinv = NULL, warn = mclust.options("warn"), 
          x = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the BIC is to be calculated. 
    The default is \code{G=1:9}, unless the argument \code{x} is specified, 
    in which case the default is taken from the values associated 
    with \code{x}. 
  }
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. The help file for
    \code{\link{mclustModelNames}} describes the available models.
    The default is:
    \describe{
        \item{\code{c("E", "V")}}{for univariate data}
        \item{\code{mclust.options("emModelNames")}}{for multivariate data (n > d)}
        \item{\code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}}{the spherical and diagonal models for multivariate data (n <= d)}
     }
    unless the argument \code{x} is specified, in which case
    the default is taken from the values associated with \code{x}. 
   }
  \item{prior}{
    The default assumes no prior, but this argument allows specification of a 
    conjugate prior on the means and variances through the function 
    \code{priorControl}.
    }
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
  \item{initialization}{
   A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. For multivariate data, the default is to compute 
    a hierarchical clustering tree by applying function \code{hc} with
    \code{modelName = "VVV"} to the data or a
    subset as indicated by the \code{subset} argument. 
    The hierarchical clustering results are to start EM.  
    For univariate data, the default is to use quantiles to start EM.  
  }
  \item{\code{subset}}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.
  }
  \item{\code{noise}}{
    A logical or numeric vector indicating an initial guess as to
    which observations are noise in the data. If numeric the entries
    should correspond to row indexes of the data. If supplied, a noise
    term will be added to the model in the estimation.
  }
  }
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region.
    The default is determined by applying function \code{hypvol} to the data. 
    Used only if an initial guess as to which observations are noise 
    is supplied.
  }
  \item{warn}{
     A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when
     estimation fails. 
     The default is controlled by \code{\link{mclust.options}}.
  }
  \item{x}{
     An object of class \code{'mclustBIC'}. If supplied, \code{mclustBIC}
     will use the settings in \code{x} to produce another object of
     class \code{'mclustBIC'}, but with \code{G} and \code{modelNames}
     as specified in the arguments. Models that have already been computed
     in \code{x} are not recomputed. All arguments to \code{mclustBIC} 
     except \code{data}, \code{G} and \code{modelName} are
     ignored and their values are set as specified in the attributes of
     \code{x}. 
     Defaults for \code{G} and \code{modelNames} are taken from \code{x}.
  }
  \item{\dots }{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}

\value{
Return an object of class \code{'mclustBIC'} containing the Bayesian Information
Criterion for the specified mixture models numbers of clusters. 
Auxiliary information returned as attributes.

The corresponding \code{print} method shows the matrix of values and the top models according to the BIC criterion.
}

\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  
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
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustModel}}, 
  \code{\link{summary.mclustBIC}}, 
  \code{\link{hc}},
  \code{\link{me}},
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
irisBIC
plot(irisBIC)

\dontrun{
subset <- sample(1:nrow(iris), 100)
irisBIC <- mclustBIC(iris[,-5], initialization=list(subset = subset))
irisBIC
plot(irisBIC)

irisBIC1 <- mclustBIC(iris[,-5], G=seq(from=1,to=9,by=2), 
                    modelNames=c("EII", "EEI", "EEE"))
irisBIC1
plot(irisBIC1)
irisBIC2  <- mclustBIC(iris[,-5], G=seq(from=2,to=8,by=2), 
                       modelNames=c("VII", "VVI", "VVV"), x= irisBIC1)
irisBIC2
plot(irisBIC2)
}

nNoise <- 450
set.seed(0)
poissonNoise <- apply(apply( iris[,-5], 2, range), 2, function(x, n) 
                      runif(n, min = x[1]-.1, max = x[2]+.1), n = nNoise)
set.seed(0)
noiseInit <- sample(c(TRUE,FALSE),size=nrow(iris)+nNoise,replace=TRUE,
                    prob=c(3,1))
irisNdata <- rbind(iris[,-5], poissonNoise)
irisNbic <- mclustBIC(data = irisNdata, G = 1:5,
                      initialization = list(noise = noiseInit))
irisNbic
plot(irisNbic)
}
\keyword{cluster}
% docclass is function
