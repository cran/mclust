\name{MclustBootstrap}
\alias{MclustBootstrap}
\alias{print.MclustBootstrap}

\title{Bootstrap Inference for Gaussian finite mixture models}

\description{Bootstrap estimation of standard errors and percentile bootstrap confidence intervals for the parameters of a Gaussian mixture model.}

\usage{
MclustBootstrap(object, nboot = 999, type = c("nonpara", "wlb"), 
                verbose = TRUE, \dots)
}

\arguments{
  \item{object}{An object of class \code{'Mclust'} or \code{'densityMclust'} providing an estimated Gaussian mixture model.}
  \item{nboot}{The number of bootstrap replications.}
  \item{type}{A character string specifying the type of bootstrap to use: 
  \tabular{ll}{
    \code{"nonpara"} \tab = nonparametric bootstrap \cr
    \code{"wlb"} \tab = weighted likelihood bootstrap
    }
  }
  \item{verbose}{A logical, if \code{TRUE} and the session is interactive a text progress bar is displayed during the bootstrap procedure.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{Nonparametric bootstrap is used to compute the bootstrap distribution of mixture parameters. For a fitted mixture model given by the first argument of the function call with \code{object$G} mixture components and covariances parameterisation model \code{object$modelName}, the bootstrap procedure generates \code{nboot} bootstrap samples of the same size as the original data by re-sampling with replacement from the observed data.
The resulting bootstrap distribution can then be used to obtain standard errors and percentile confidence intervals by the use of \code{\link{summary.MclustBootstrap}} function.}

\value{An object of class \code{'MclustBootstrap'} with the following components:

  \item{G}{A value specifying the number of mixture components.} 
  \item{modelName}{A character string specifying the mixture model covariances 
    parameterisation (see \code{\link{mclustModelNames}}).}
  \item{parameters}{A list of estimated parameters for the mixture components with the following components:  
  \describe{
     \item{\code{pro}}{a vector of mixing proportions.}
     \item{\code{mean}}{a matrix of means for each component.}
     \item{\code{variance}}{an array of covariance matrices for each component.}
     }
  }
  \item{nboot}{The number of bootstrap replications.}
  \item{type}{The type of bootstrap computation employed.}
  \item{pro}{A matrix of dimension (\code{nboot} x \code{G}) containing the 
    bootstrap distribution for the mixing proportion.}
  \item{mean}{An array of dimension (\code{nboot} x \code{d} x \code{G}), 
    where \code{d} is the dimension of the data, containing the bootstrap 
    distribution for the component means.}
  \item{variance}{An array of dimension (\code{nboot} x \code{d} x \code{d} x 
    \code{G}), where \code{d} is the dimension of the data, containing the 
    bootstrap distribution for the component covariances.}
}

\references{
Davison, A. and Hinkley, D. (1997) \emph{Bootstrap Methods and Their Applications}. Cambridge University Press.

McLachlan, G.J. and Peel, D. (2000) \emph{Finite Mixture Models}. Wiley.
}

\seealso{\code{\link{summary.MclustBootstrap}}, \code{\link{plot.MclustBootstrap}}, \code{\link{Mclust}}, \code{\link{densityMclust}}.}

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
