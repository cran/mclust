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
\alias{simEVE}
\alias{simEVV}
\alias{simVEE}
\alias{simVVE}

\title{
  Simulate from a Parameterized MVN Mixture Model
}
\description{
  Simulate data from a parameterized MVN mixture model.
}
\usage{
simE(parameters, n, seed = NULL, \dots)
simV(parameters, n, seed = NULL, \dots)
simEII(parameters, n, seed = NULL, \dots)
simVII(parameters, n, seed = NULL, \dots)
simEEI(parameters, n, seed = NULL, \dots)
simVEI(parameters, n, seed = NULL, \dots)
simEVI(parameters, n, seed = NULL, \dots)
simVVI(parameters, n, seed = NULL, \dots)
simEEE(parameters, n, seed = NULL, \dots)
simEEV(parameters, n, seed = NULL, \dots)
simVEV(parameters, n, seed = NULL, \dots)
simVVV(parameters, n, seed = NULL, \dots)
simEVE(parameters, n, seed = NULL, \dots)
simEVV(parameters, n, seed = NULL, \dots)
simVEE(parameters, n, seed = NULL, \dots)
simVVE(parameters, n, seed = NULL, \dots)
}
\arguments{
   \item{parameters}{
     A list with the following components:
     \describe{
        \item{\code{pro}}{
              A vector whose \emph{k}th component is the mixing proportion for
              the \emph{k}th component of the mixture model.
              If missing, equal proportions are assumed.
        }
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}}
              for details.
        }
      }
  }
  \item{n}{
      An integer specifying the number of data points to be simulated.
  }
  \item{seed}{
    An optional integer argument to \code{set.seed} for reproducible
    random class assignment. By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{sim}.
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
      \code{"modelName"} A character string indicating the variance 
      model used for the simulation.
  }
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep}, \code{em}
  \code{me}, \code{Mclust}, to be passed directly without the need
  to specify individual parameters as arguments. 
}

\seealso{
  \code{\link{sim}},
  \code{\link{Mclust}},
  \code{\link{mstepE}},
  \code{\link{mclustVariance}}.
}
\examples{
\dontrun{
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
simdat <- simEEV( n = 200, 
                  parameters = list(pro=c(1,1),mean=mu,variance=variance),
                  seed = NULL)

cl <- simdat[,1]

sigma <- array(apply(O, 3, function(x,y) crossprod(x*y), 
                 y = sqrt(scale*shape)), c(2,2,2))
paramList <- list(mu = mu, sigma = sigma)
coordProj( simdat, paramList = paramList, classification = cl)
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
