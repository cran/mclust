\name{dens}
\alias{dens}
\title{
  Density for Parameterized MVN Mixtures
}
\description{
  Computes densities of obseravations in parameterized MVN mixtures.
}
\usage{
dens(modelName, data, mu, logarithm, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. Possible models: \cr\cr
    "E" for spherical, equal variance (one-dimensional) \cr
    "V" for spherical, variable variance (one-dimensional) \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume, equal shape \cr
    "VEI": diagonal, varying volume, equal shape \cr
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume, varying shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation\cr 
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation \cr\cr
    For fitting a single Gaussian,\cr\cr
    "X": one-dimensional \cr
    "XII": spherical \cr
    "XXI": diagonal \cr
    "XXX": ellipsoidal
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{mu}{
    The mean for each component. If there is more than one component,
    \code{mu} is a matrix whose columns are the means of the components.
  }
  \item{logarithm}{Return logarithm of the density, rather than the
    density itself. Default: FALSE}
  \item{\dots}{
    Other arguments, such as an argument describing the variance.
    See \code{\link{cdens}}.
  }
}
\value{
  A numeric vector whose \emph{i}th component is the density of observation 
  \emph{i} in the MVN mixture specified by \code{mu} and
  \code{\dots}.
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
  \code{\link{grid1}},
  \code{\link{cdens}},
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

set.seed(0)
I <- sample(1:(2*n))
x <- x[I, ]
xclass <- xclass[I]

odd <- seq(1, 2*n, by = 2)
oddBic <- EMclust(x[odd, ]) 
oddSumry <- summary(oddBic, x[odd, ]) ## best parameter estimates
names(oddSumry)

oddDens <- dens(modelName = oddSumry$modelName, data = x,
   mu = oddSumry$mu, decomp = oddSumry$decomp, pro = oddSumry$pro)
oddDens <- do.call("dens", c(list(data = x), oddSumry))  ## alternative call

even <- odd + 1
evenBic <- EMclust(x[even, ]) 
evenSumry <- summary(evenBic, x[even, ]) ## best parameter estimates
evenDens <- do.call( "dens", c(list(data = x), evenSumry))

cbind(class = xclass, odd = oddDens, even = evenDens)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
