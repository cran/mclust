\name{summary.MclustBootstrap}
\alias{summary.MclustBootstrap}
\alias{print.summary.MclustBootstrap}

\title{Summary Function for Bootstrap Inference for Gaussian Finite Mixture Models}

\description{Summary of bootstrap distribution for the parameters of a Gaussian mixture model providing either standard errors or percentile bootstrap confidence intervals.}

\usage{
\method{summary}{MclustBootstrap}(object, what = c("se", "ci", "ave"), conf.level = 0.95, \dots)
}

\arguments{
  \item{object}{An object of class \code{'MclustBootstrap'} as returned by \code{\link{MclustBootstrap}}.}
  \item{what}{A character string: \code{"se"} for the standard errors; \code{"ci"} for the confidence intervals; \code{"ave"} for the averages.}
  \item{conf.level}{A value specifying the confidence level of the interval.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{For details about the procedure used to obtain the bootstrap distribution see \code{\link{MclustBootstrap}}.}

%\value{}

\seealso{\code{\link{MclustBootstrap}}.}

\examples{
\donttest{
data(diabetes)
X = diabetes[,-1]
modClust = Mclust(X) 
bootClust = MclustBootstrap(modClust)
summary(bootClust, what = "se")
summary(bootClust, what = "ci")

data(acidity)
modDens = densityMclust(acidity)
modDens = MclustBootstrap(modDens)
summary(modDens, what = "se")
summary(modDens, what = "ci")
}
}

\keyword{htest}
\keyword{cluster}
