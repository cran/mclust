\name{mclustDAtrain}
\alias{mclustDAtrain}
\title{
  MclustDA Training
}
\description{
  Training phase for MclustDA discriminant analysis.
}
\usage{
mclustDAtrain(data, labels, G, emModelNames, eps, tol, itmax, 
              equalPro, warnSingular, verbose)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{labels}{
    A numeric or character vector assigning a class label to each observation.
  }
  \item{G}{
    An integer vector specifying the numbers of Gaussian mixture
    components (clusters) for which the BIC is to be calculated (the
    same specification is used for all classes). Default: \emph{1:9}.
  }
  \item{emModelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. Possible models: \cr
    "E" for spherical, equal variance (one-dimensional)\cr 
    "V" for spherical, variable variance (one-dimensional)\cr 
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
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due
    to computational singularity in covariances. Smaller values of
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
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is \code{warnSingular=FALSE}.
  }
  \item{verbose}{
    A logical value indicating whether or not to print the models and
    numbers of components for each class.
    Default:\code{verbose=TRUE}.
  }
}
\value{
  A list in which each element gives the optimal parameters for the model
  best fitting each class according to BIC.
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
  \code{\link{summary.mclustDAtrain}}, 
  \code{\link{mclustDAtest}}, 
  \code{\link{EMclust}}, 
  \code{\link{hc}},
  \code{\link{mclustOptions}}
}
\examples{
n <- 250 ## create artificial data
set.seed(0)
par(pty = "s")
x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))
\dontrun{
mclust2Dplot(x, classification = xclass, type="classification", ask=FALSE)}

odd <- seq(1, 2*n, 2)
train <- mclustDAtrain(x[odd, ], labels = xclass[odd]) ## training step
summary(train)

even <- odd + 1
test <- mclustDAtest(x[even, ], train) ## compute model densities
clEven <- summary(test)$class ## classify training set
compareClass(clEven,xclass[even])
}
\keyword{cluster}
% docclass is function
