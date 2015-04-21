\name{mclustModelNames}
\alias{mclustModelNames}
\title{
   MCLUST Model Names 
}
\description{
Description of model names used in the \emph{MCLUST} package.
}
\usage{
mclustModelNames(model)
}
\arguments{
  \item{model}{A string specifying the model.}
}

\details{
The following models are available in package \pkg{mclust}:
\tabular{lcl}{
\bold{univariate mixture} \tab \tab \cr
    \code{"E"} \tab = \tab equal variance  (one-dimensional) \cr
    \code{"V"} \tab = \tab variable variance (one-dimensional) \cr
\bold{multivariate mixture}  \tab \tab \cr
    \code{"EII"} \tab = \tab spherical, equal volume \cr
    \code{"VII"} \tab = \tab spherical, unequal volume \cr
    \code{"EEI"} \tab = \tab diagonal, equal volume and shape\cr 
    \code{"VEI"} \tab = \tab diagonal, varying volume, equal shape\cr 
    \code{"EVI"} \tab = \tab diagonal, equal volume, varying shape \cr
    \code{"VVI"} \tab = \tab diagonal, varying volume and shape \cr
    \code{"EEE"} \tab = \tab ellipsoidal, equal volume, shape, and orientation \cr
    \code{"EVE"} \tab = \tab ellipsoidal, equal volume and orientation (*)\cr
    \code{"VEE"} \tab = \tab ellipsoidal, equal shape and orientation (*)\cr
    \code{"VVE"} \tab = \tab ellipsoidal, equal orientation (*)\cr
    \code{"EEV"} \tab = \tab ellipsoidal, equal volume and equal shape\cr
    \code{"VEV"} \tab = \tab ellipsoidal, equal shape \cr
    \code{"EVV"} \tab = \tab ellipsoidal, equal volume (*)\cr
    \code{"VVV"} \tab = \tab ellipsoidal, varying volume, shape, and orientation \cr
\bold{single component}  \tab \tab \cr
    \code{"X"} \tab = \tab univariate normal \cr
    \code{"XII"} \tab = \tab spherical multivariate normal \cr
    \code{"XXI"} \tab = \tab diagonal multivariate normal \cr
    \code{"XXX"} \tab = \tab ellipsoidal multivariate normal 
}
(*) new models in \pkg{mclust} version >= 5.0.0.
}

\value{Returns a list with the following components:
\item{model}{a character string indicating the model (as in input).}
\item{type}{the description of the indicated model (see Details section).}
}

\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}.

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{Mclust}},
  \code{\link{mclustBIC}}
}
\examples{
mclustModelNames("E")
mclustModelNames("EEE")
mclustModelNames("VVV")
mclustModelNames("XXI")
}
\keyword{cluster}

