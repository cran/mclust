\name{MclustBootstrap}
\alias{MclustBootstrap}
\alias{print.MclustBootstrap}

\title{Bootstrap Inference for Gaussian finite mixture models}

\description{Boostrap estimation of standard errors and percentile boostrap confidence intervals for the parameters of a Gaussian mixture model.}

\usage{
MclustBootstrap(object, nboot = 999, verbose = TRUE, \dots)
}

\arguments{
  \item{object}{An object of class \code{'Mclust'} or \code{'densityMclust'} providing an estimated Gaussian mixture model.}
  \item{nboot}{The number of bootstrap replications.}
  \item{verbose}{A logical, if \code{TRUE} and the session is interactive a text progress bar is displayed during the bootstrap procedure.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{Nonparametric boostrap is used to compute the boostrap distribution of mixture parameters. For a fitted mixture model given by the first argument of the function call with \code{object$G} mixture components and covariances parameterisation model \code{object$modelName}, the boostrap procedure generates \code{nboot} bootstrap samples of the same size as the original data by re-sampling with replacement from the observed data.
The resulting bootstrap distribution can then be used to obtain standard errors and percentile confidence intervals by the use of \code{\link{summary.MclustBootstrap}} function.}

\value{An object of class \code{'MclustBootstrap'} with the following components:

  \item{G}{A value specifying the number of mixture components.} 
  \item{modelName}{A character string specifying the mixture model covariances 
    parameterisation (see \code{\link{mclustModelNames}}).}
  \item{nboot}{The number of bootstrap replications.}
  \item{pro}{A matrix of dimension \code{nboot} x \code{G} containing the 
    boostrap distribution for the mixing proportion.}
  \item{mean}{An array of dimension \code{nboot} x \code{d} x \code{G}, 
    where \code{d} is the dimension of the data, containing the boostrap 
    distribution for the component means.}
  \item{variance}{An array of dimension \code{nboot} x \code{d} x \code{d} x \code{G}, 
    where \code{d} is the dimension of the data, containing the boostrap 
    distribution for the component covariances.}
}

\references{
Davison, A. and Hinkley, D. (1997) \emph{Bootstrap Methods and Their Applications}. Cambridge University Press.

McLachlan, G.J. and Peel, D. (2000) \emph{Finite Mixture Models}. Wiley.
}

\seealso{\code{\link{Mclust}}, \code{\link{densityMclust}}.}

\examples{
\dontrun{
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
