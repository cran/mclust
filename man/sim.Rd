\name{sim}
\alias{sim}
\title{
  Simulate from Parameterized MVN Mixture Models
}
\description{
  Simulate data from parameterized MVN mixture models.
}
\usage{
sim(modelName, mu, \dots, seed = 0)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. Possible models: \cr\cr
    "E": equal variance  (one-dimensional) \cr
    "V": variable variance (one-dimensional) \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume, equal shape \cr
    "VEI": diagonal, varying volume, equal shape \cr
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume, varying shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
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
	  ellipsoidal models ("EEV", "VEV"). This is a list described in
	  \code{\link{cdens}}. 
	}
	\item{Sigma}{
	  for the equal variance model "EEE". A \emph{d} by \emph{d}
	  matrix giving the common covariance for all
	  components of the mixture model.
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
    }
    \item{pro}{
      Component mixing proportions. If missing, equal proportions are
      assumed. 
    }
    \item{n}{
      An integer specifying the number of data points to be simulated.
    }
  }
  \item{seed}{
    A integer between 0 and 1000, inclusive, for specifying a seed for 
    random class assignment. The default value is 0.
  }
}
\value{
  A data set consisting of n points simulated from the specified 
  MVN mixture model.
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep}, \code{em},
  \code{me}, or \code{EMclust} to be passed directly without the need to
  specify individual parameters as arguments. 
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
  \code{\link{simE}}, \ldots,
  \code{\link{simVVV}},
  \code{\link{EMclust}},
  \code{\link{mstep}},
  \code{\link{do.call}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])

irisBic <- EMclust(irisMatrix)
irisSumry <- summary(irisBic,irisMatrix)
names(irisSumry)
irisSim <- sim(modelName = irisSumry$modelName, n = dim(irisMatrix)[1],
       mu = irisSumry$mu, decomp = irisSumry$decomp, pro = irisSumry$pro)
\dontrun{
irisSim <- do.call("sim", irisSumry) ## alternative call
}

par(pty = "s", mfrow = c(1,2))
dimens <- c(1,2)
xlim <- range(rbind(irisMatrix,irisSim)[,dimens][,1])
ylim <- range(rbind(irisMatrix,irisSim)[,dimens][,2])

cl <- irisSumry$classification
coordProj(irisMatrix, par=irisSumry, classification=cl, dimens=dimens, 
          xlim=xlim, ylim=ylim)
cl <- attr(irisSim,"classification")
coordProj(irisSim, par=irisSumry, classification=cl, dimens=dimens, 
          xlim=xlim, ylim=ylim)

irisSumry3 <- summary(irisBic,irisMatrix, G=3)
irisSim3 <- do.call("sim", c(list(n = 500, seed = 1), irisSumry3))
clPairs(irisSim3, cl = attr(irisSim3,"classification"))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
