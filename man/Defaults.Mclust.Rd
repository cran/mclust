\name{Defaults.Mclust}
\alias{.Mclust}
\title{
  List of values controlling defaults for some MCLUST functions.
}
\description{
  A named list of values including tolerances for singularity and
  convergence assessment, and an enumeration of models used as defaults
  in MCLUST functions. 
}
\value{
  A list with the following components:
  
  \item{eps}{
    A scalar tolerance for deciding when to terminate computations due
    to computational singularity in covariances. Smaller values of
    \code{eps} allow computations to proceed nearer to singularity. The
    default is the relative machine precision
    \code{.Machine\$double.eps}, which is approximately \$2e-16\$ on
    IEEE-compliant machines. 
  }
  \item{tol}{
    A vector of length two giving relative convergence tolerances for the 
    loglikelihood and for parameter convergence in the inner loop for models
    with iterative M-step ("VEI", "VEE", "VVE", "VEV"), respectively.
    The default is \code{c(1.e-5,1.e-5)}.
  }
  \item{itmax}{
    A vector of length two giving integer limits on the number of EM
    iterations and on the number of iterations in the inner loop for
    models with iterative M-step ("VEI", "VEE", "VVE", "VEV"),
    respectively. The default is \code{c(Inf,Inf)} allowing termination
    to be completely governed by \code{tol}.
  }
  \item{equalPro}{
    Logical variable indicating whether or not the mixing proportions are
    equal in the model. Default: \code{equalPro = FALSE}.
  }
  \item{warnSingular}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered. Default: \code{warnSingular = TRUE}.
  }
  \item{emModelNames}{
    A vector of character strings indicating the models to be used for
    multivariate data in the functions such as \code{EMclust}
    and \code{mclustDAtrain} that involve multiple models. The default
    is all of the multivariate models available in MCLUST: \cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape \cr
    "VEI": diagonal, varying volume, equal shape \cr
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
  }
  \item{hcModelName}{
    A vector of two character strings 
    giving the name of the model to be used in the hierarchical
    clustering phase for univariate and multivariate data, respectively,
    in \code{EMclust} and \code{EMclustN}. The default is
    \code{c("V","VVV")}, giving the unconstrained model in each case. 
  }
  \item{symbols}{A vector whose entries are either integers
    corresponding to graphics symbols or single characters for plotting
    for classifications. Classes are assigned symbols in the given
    order. 
  }
}
\details{
  A function \code{mclustOptions} is supplied for assigning values to
  the \code{.Mclust} list. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association}. 
  See \url{http://www.stat.washington.edu/tech.reports} (No. 380, 2000). 
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/tech.reports}.
}
\seealso{
  \code{\link{mclustOptions}},
  \code{\link{EMclust}},
  \code{\link{mclustDAtrain}},
  \code{\link{em}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mstep}}
}
\examples{
data(iris)
irisMatrix <- iris[,1:4]
irisClass <- iris[,5]

.Mclust
.Mclust <- mclustOptions(tol = 1.e-6, emModelNames = c("VII", "VVI", "VVV"))
.Mclust
irisBic <- EMclust(irisMatrix)
summary(irisBic, irisMatrix)
.Mclust <- mclustOptions() # restore defaults
.Mclust
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
