\name{cdensE}
\alias{cdensE}
\alias{cdensV}
\alias{cdensEII}
\alias{cdensVII}
\alias{cdensEEI}
\alias{cdensVEI}
\alias{cdensEVI}
\alias{cdensVVI}
\alias{cdensEEE}
\alias{cdensEEV}
\alias{cdensVEV}
\alias{cdensVVV}
\title{
  Component Density for a Parameterized MVN Mixture Model
}
\description{
  Computes component densities for points in a parameterized MVN mixture model.
}
\usage{
cdensE(data, mu, sigmasq, eps, warnSingular, logarithm = FALSE, \dots)
cdensV(data, mu, sigmasq, eps, warnSingular, logarithm = FALSE, \dots)
cdensEII(data, mu, sigmasq, eps, warnSingular, logarithm = FALSE, \dots)
cdensVII(data, mu, sigmasq, eps, warnSingular, logarithm = FALSE, \dots)
cdensEEI(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensVEI(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensEVI(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensVVI(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensEEE(data, mu, eps, warnSingular, logarithm = FALSE, \dots)
cdensEEV(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensVEV(data, mu, decomp, eps, warnSingular, logarithm = FALSE, \dots)
cdensVVV(data, mu, eps, warnSingular, logarithm = FALSE, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the
    means of the components.
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
  \item{logarithm}{
    A logical value indicating whether or not the logarithm of the
    component densities should be returned.
    The default is to return the component densities,
    obtained from the log component densities by exponentiation.
  }
  \item{\dots}{
    An argument giving the variance that takes one of the following forms:
    \describe{
      \item{decomp}{for models "EII" and "VII"; see above.}
      \item{cholSigma}{see Sigma, for "EEE".}
      \item{Sigma}{ 
	for the equal variance model "EEE". A \emph{d} by \emph{d} matrix
	giving the common covariance for all components of the mixture
	model.
      }
      \item{cholsigma}{see sigma, for "VVV".}
      \item{sigma}{ 
	for the unconstrained variance model "VVV".
	A \emph{d} by \emph{d} by \emph{G} matrix array whose 
	\code{[,,k]}th entry is the covariance matrix for the
	\emph{k}th component of the mixture model.
      }
      
      The form of the variance specification is the same as for the output
      for the \code{em}, \code{me}, or \code{mstep} methods for the
      specified mixture model.
      
      Also used to catch unused arguments from a \code{do.call} call.
    }
  }
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due
    to computational singularity in covariances. Smaller values of
    \code{eps} allow computations to proceed nearer to singularity. The
    default is \code{.Mclust\$eps}. 
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered. The default is
    \code{.Mclust\$warnSingular}. 
  }
}
\value{
  A numeric matrix whose \code{[i,j]}th 
  entry is the density of observation \emph{i} in component \emph{j}. 
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
  \code{\link{cdens}},
  \code{\link{dens}},
  \code{\link{EMclust}},
  \code{\link{mstep}},
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

modelVII <- meVII(x, z = unmap(xclass))
modelVVI <- meVVI(x, z = unmap(xclass))
modelVVV <- meVVV(x, z = unmap(xclass))

names(modelVII)
args(cdensVII)
cdenVII <- cdensVII(data = x, mu = modelVII$mu, pro = modelVII$pro,
                    decomp = modelVII$decomp)
names(modelVVI)
args(cdensVVI)
cdenVVI <- cdensVII(data = x, mu = modelVVI$mu, pro = modelVVI$pro,
                    decomp = modelVVI$decomp)
names(modelVVV)
args(cdensVVV)
cdenVVV <- cdensVVV( data = x, mu = modelVVV$mu, pro = modelVVV$pro,
                     cholsigma = modelVVV$cholsigma)

cbind(class=xclass,VII=map(cdenVII),VVI=map(cdenVVI),VVV=map(cdenVVV))

## alternative call

\dontrun{
cdenVII <- do.call("cdensVII", c(list(data = x), modelVII))
cdenVVI <- do.call("cdensVVI", c(list(data = x), modelVVI))
cdenVVV <- do.call("cdensVVV", c(list(data = x), modelVVV))

cbind(class=xclass,VII=map(cdenVII),VVI=map(cdenVVI),VVV=map(cdenVVV))
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
