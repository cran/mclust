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
simE(parameters, n, seed = 0, \dots)
simV(parameters, n, seed = 0, \dots)
simEII(parameters, n, seed = 0, \dots)
simVII(parameters, n, seed = 0, \dots)
simEEI(parameters, n, seed = 0, \dots)
simVEI(parameters, n, seed = 0, \dots)
simEVI(parameters, n, seed = 0, \dots)
simVVI(parameters, n, seed = 0, \dots)
simEEE(parameters, n, seed = 0, \dots)
simEEV(parameters, n, seed = 0, \dots)
simVEV(parameters, n, seed = 0, \dots)
simVVV(parameters, n, seed = 0, \dots)
}
\arguments{
   \item{parameters}{
     A list with the following components:
     \describe{
        \item{pro}{
              A vector whose \emph{k}th component is the mixing proportion for
              the \emph{k}th component of the mixture model.
              If missing, equal proportions are assumed.
        }
        \item{mean}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{variance}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{mclustVariance}
              for details.
        }
      }
  }
  \item{n}{
      An integer specifying the number of data points to be simulated.
  }
  \item{seed}{
    A integer for specifying a seed for random class assignment. 
    The default value is 0.
  }
  \item{\dots }{
     Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A matrix in which first column is the classification and the remaining
  columns are the \code{n} observations simulated from the specified MVN 
  mixture model.
\item{Attributes:}{
    \itemize{
      \item \code{"modelName"}
       A character string indicating the variance model used for
       the simulation.
    }
  }
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep}, \code{em}
  \code{me}, \code{Mclust}, to be passed directly without the need
  to specify individual parameters as arguments. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{sim}},
  \code{\link{Mclust}},
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

variance <- list(d= d, G = G, scale = scale, shape = shape, orientation = O)
mu <- matrix(0, d, G) ## center at the origin
simdat <- simEEV(n=200,parameters=list(pro=c(1,1),mean=mu,variance=variance))

cl <- simdat[,1]

\dontrun{
sigma <- array(apply(O, 3, function(x,y) crossprod(x*y), 
                 y = sqrt(scale*shape)), c(2,2,2))
paramList <- list(mu = mu, sigma = sigma)
coordProj( simdat, paramList = paramList, classification = cl)
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
