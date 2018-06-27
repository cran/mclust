\name{MclustBootstrap}
\alias{MclustBootstrap}
\alias{print.MclustBootstrap}

\title{Resampling-based Inference for Gaussian finite mixture models}

\description{Bootstrap or jackknife estimation of standard errors and percentile bootstrap confidence intervals for the parameters of a Gaussian mixture model.}

\usage{
MclustBootstrap(object, nboot = 999, type = c("bs", "wlbs", "pb", "jk"),
                max.nonfit = 10*nboot, verbose = interactive(), 
                \dots)
}

\arguments{
  \item{object}{An object of class \code{'Mclust'} or \code{'densityMclust'} providing an estimated Gaussian mixture model.}
  \item{nboot}{The number of bootstrap replications.}
  \item{type}{A character string specifying the type of resampling to use: 
  \tabular{ll}{
    \code{"bs"} \tab = nonparametric bootstrap \cr
    \code{"wlbs"} \tab = weighted likelihood bootstrap \cr
    \code{"pb"} \tab = parametric bootstrap \cr
    \code{"jk"} \tab = jackknife \cr
    }
  }
  \item{max.nonfit}{The maximum number of non-estimable models allowed.}
  \item{verbose}{A logical controlling if a text progress bar is displayed during the resampling procedure. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{For a fitted Gaussian mixture model with \code{object$G} mixture components and covariances parameterisation \code{object$modelName}, this function returns either the boostrap distribution or the jackknife distribution of mixture parameters. In the former case, the nonparametric bootstrap or the weighted likelihood bootstrap approach could be used, so the the bootstrap procedure generates \code{nboot} bootstrap samples of the same size as the original data by resampling with replacement from the observed data. In the jackknife case, the procedure considers all the samples obtained by omitting one observation at time.

The resulting resampling distribution can then be used to obtain standard errors and percentile confidence intervals by the use of \code{\link{summary.MclustBootstrap}} function.}

\value{An object of class \code{'MclustBootstrap'} with the following components:

  \item{n}{The number of observations in the data.}
  \item{d}{The dimension of the data.}
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
  \item{nboot}{The number of bootstrap replications if \code{type = "bs"} or \code{type = "wlbs"}. The sample size if \code{type = "jk"}.}
  \item{type}{The type of resampling approach used.}
  \item{nonfit}{The number of resamples that did not convergence during the procedure.}
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

O'Hagan A., Murphy T. B., Gormley I. C. and Scrucca L. (2015) On Estimation of Parameter Uncertainty in Model-Based Clustering. Submitted to \emph{Computational Statistics}. 
}

\seealso{\code{\link{summary.MclustBootstrap}}, \code{\link{plot.MclustBootstrap}}, \code{\link{Mclust}}, \code{\link{densityMclust}}.}

\examples{
\dontrun{
data(diabetes)
X <- diabetes[,-1]
modClust <- Mclust(X) 
bootClust <- MclustBootstrap(modClust)
summary(bootClust, what = "se")
summary(bootClust, what = "ci")

data(acidity)
modDens <- densityMclust(acidity)
modDens <- MclustBootstrap(modDens)
summary(modDens, what = "se")
summary(modDens, what = "ci")
}
}

\keyword{htest}
\keyword{cluster}
