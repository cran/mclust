\name{mclustDAtrainN}
\alias{mclustDAtrainN}
\title{MclustDA training with noise}
\description{
  Training phase for MclustDA discriminant analysis in the presence
  of noise.
}
\usage{
mclustDAtrainN(data, labels, emModelNames, G, hcModelName,
               equalPro = FALSE, noise, Vinv)
}
\arguments{
  \item{data}{A numeric vector, matrix, or data frame of
    observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables.
  }
  \item{labels}{A numeric or character vector assigning a class label to
    each observation.
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
  \item{G}{
    An integer vector specifying the numbers of Gaussian mixture
    components (clusters) for which the BIC is to be calculated (the
    same specification is used for all classes). Default: \emph{0:9}.
  }
  \item{hcModelName}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. The default is to compute a hierarchical
    clustering tree by applying function \code{hc} with
    \code{modelName = .Mclust\$hcModelName[1]} to univariate data and
    \code{modelName = .Mclust\$hcModelName[2]} to multivariate data or a
    subset as indicated by the \code{subset} argument. The hierarchical
    clustering results are used as starting values for EM. 
  }
  \item{equalPro}{
    Logical variable indicating whether or not the mixing
    proportions are equal in the model. The default is
    \code{.Mclust\$equalPro}.
  }
  \item{noise}{A logical vector indicating whether or not observations
    are initially estimated to noise in the data. If there is no noise
    \code{mclustDAtrain} should be use rather than \code{mclustDAtrainN}. 
  }
  \item{Vinv}{
    An estimate of the reciprocal hypervolume of the data region.
    The default is determined by applying function
    \code{hypvol} to the data.
  }
}
\value{
  An object of class \code{mclustDAtrain}. See \code{\link{mclustDAtrain}}.
}
\references{  C. Fraley and A. E. Raftery (2002a).
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
  \code{\link{mclustDAtrain}},
  \code{\link{mclustDA}},
  \code{\link{mclustDAtest}}
}
\keyword{cluster}

