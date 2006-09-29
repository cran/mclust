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
cdensE(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensV(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensEII(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensVII(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensEEI(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensVEI(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensEVI(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensVVI(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensEEE(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensEEV(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensVEV(data, logarithm = FALSE, parameters, warn = NULL, \dots)
cdensVVV(data, logarithm = FALSE, parameters, warn = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{logarithm}{
    A logical value indicating whether or not the logarithm of the
    component densities should be returned.
    The default is to return the component densities,
    obtained from the log component densities by exponentiation.
  }
 \item{parameters}{
     The parameters of the model:
     \describe{
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
 \item{warn}{
    A logical value indicating whether or not a warning should be issued
    when computations fail. The default is \code{warn=FALSE}.
  }
 \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A numeric matrix whose \code{[i,j]}th 
  entry is the density of observation \emph{i} in component \emph{j}. 
  The densities are not scaled by mixing proportions.
}
\note{
  When one or more component densities are very large in magnitude,
  then it may be possible to compute the logarithm of the component
  densities but not the component densities themselves due to overflow.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{cdens}},
  \code{\link{dens}},
  \code{\link{mclustBIC}},
  \code{\link{mstep}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{

z2 <- unmap(hclass(hcVVV(faithful),2)) # initial value for 2 class case

model <- meVVV(data=faithful, z=z2)
cdensVVV(data=faithful, logarithm = TRUE, parameters = model$parameters)

z2 <- unmap(cross[,1])

model <- meEEV(data = cross[,-1], z = z2)

EEVdensities <- cdensEEV( data = cross[,-1], parameters = model$parameters)

cbind(cross[,-1],map(EEVdensities))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
