\name{simE}
\alias{simE}
\alias{simV}
\alias{simEII}
\alias{simVII}
\alias{simEEI}
\alias{simVEI}
\alias{simEVI}
\alias{simVVI}
\alias{simEEE}
\alias{simEEV}
\alias{simVEV}
\alias{simVVV}
\title{
  Simulate from a Parameterized MVN Mixture Model
}
\description{
  Simulate data from a parameterized MVN mixture model.
}
\usage{
simE(mu, sigmasq, pro, \dots, seed = 0)
simV(mu, sigmasq, pro, \dots, seed = 0)
simEII(mu, sigmasq, pro, \dots, seed = 0)
simVII(mu, sigmasq, pro, \dots, seed = 0)
simEEI(mu, decomp, pro, \dots, seed = 0)
simVEI(mu, decomp, pro, \dots, seed = 0)
simEVI(mu, decomp, pro, \dots, seed = 0)
simVVI(mu, decomp, pro, \dots, seed = 0)
simEEE(mu, pro, \dots, seed = 0)
simEEV(mu, decomp, pro, \dots, seed = 0)
simVEV(mu, decomp, pro, \dots, seed = 0)
simVVV mu, pro, \dots, seed = 0)
}
\arguments{
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the
    means of the components.
  }
  \item{sigmasq}{ 
    for the one-dimensional models ("E", "V") and spherical models
    ("EII", "VII"). This is either a vector whose \emph{k}th component
    is the variance for the \emph{k}th component in the mixture model
    ("V" and "VII"), or a scalar giving the common variance for all
    components in the mixture model ("E" and "EII"). 
  }
  \item{decomp}{ 
    for the diagonal models ("EEI", "VEI", "EVI", "VVI") and some
    ellipsoidal models ("EEV", "VEV"). This is a list described in
    \code{\link{cdens}}.
  }
  \item{pro}{
    Component mixing proportions. If missing, equal proportions are assumed.
  }
  \item{\dots}{
    \item Other terms describing variance:
    \describe{
      \item{Sigma}{ 
	for the equal variance model "EEE".
	A \emph{d} by \emph{d} matrix giving the common covariance for
	all components of the mixture model.
      }
      \item{sigma}{ 
	for the unconstrained variance model "VVV".
	A \emph{d} by \emph{d} by \emph{G} matrix array whose 
	\code{[,,k]}th entry is the covariance matrix for the
	\emph{k}th component of the mixture model.
      }
      The form of the variance specification is the same as for the output
      for the \code{em}, \code{me}, or \code{mstep} methods for the
      specified mixture model.  
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
  A data set consisting of \code{n} points simulated from 
  the specified MVN mixture model.
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep}, \code{em}
  \code{me}, or \code{EMclust}, to be passed directly without the need
  to specify individual parameters as arguments. 
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
  \code{\link{sim}},
  \code{\link{EMclust}},
  \code{\link{mstepE}},
  \code{\link{do.call}}
}
\examples{
d <- 2
G <- 2
scale <- 1
shape <- c(1, 9)

O1 <- diag(2)
O2 <- diag(2)[,c(2,1)]
O <- array(cbind(O1,O2), c(2, 2, 2))
O

decomp <- list(d= d, G = G, scale = scale, shape = shape, orientation = O)
mu <- matrix(0, d, G) ## center at the origin
simdat <- simEEV(n=200, mu=mu, decomp=decomp, pro = c(1,1))

cl <- attr(simdat, "classification")
sigma <- array(apply(O, 3, function(x,y) crossprod(x*y), 
                 y = sqrt(scale*shape)), c(2,2,2))
paramList <- list(mu = mu, sigma = sigma)
coordProj( simdat, paramList = paramList, classification = cl)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
