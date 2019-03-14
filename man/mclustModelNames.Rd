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
The following models are available in package \pkg{mclust}:\cr

\bold{univariate mixture} \cr
\describe{
\item{\code{"E"}}{equal variance (one-dimensional)}
\item{\code{"V"}}{variable/unqual variance (one-dimensional)}
}
\bold{multivariate mixture}\cr
\describe{
\item{\code{"EII"}}{spherical, equal volume}
\item{\code{"VII"}}{spherical, unequal volume}
\item{\code{"EEI"}}{diagonal, equal volume and shape}
\item{\code{"VEI"}}{diagonal, varying volume, equal shape}
\item{\code{"EVI"}}{diagonal, equal volume, varying shape}
\item{\code{"VVI"}}{diagonal, varying volume and shape}
\item{\code{"EEE"}}{ellipsoidal, equal volume, shape, and orientation}
\item{\code{"EVE"}}{ellipsoidal, equal volume and orientation (*)}
\item{\code{"VEE"}}{ellipsoidal, equal shape and orientation (*)}
\item{\code{"VVE"}}{ellipsoidal, equal orientation (*)}
\item{\code{"EEV"}}{ellipsoidal, equal volume and equal shape}
\item{\code{"VEV"}}{ellipsoidal, equal shape}
\item{\code{"EVV"}}{ellipsoidal, equal volume (*)}
\item{\code{"VVV"}}{ellipsoidal, varying volume, shape, and orientation}
}
\bold{single component}\cr
\describe{
\item{\code{"X"}}{univariate normal}
\item{\code{"XII"}}{spherical multivariate normal}
\item{\code{"XXI"}}{diagonal multivariate normal}
\item{\code{"XXX"}}{ellipsoidal multivariate normal}
}
(*) new models in \pkg{mclust} version >= 5.0.0.
}

\value{Returns a list with the following components:
\item{model}{a character string indicating the model (as in input).}
\item{type}{the description of the indicated model (see Details section).}
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

