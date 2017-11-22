\name{logLik.Mclust}
\alias{logLik.Mclust}

\title{Log-Likelihood of a \code{Mclust} object}

\description{
Returns the log-likelihood for a \code{'Mclust'} object.}

\usage{
  \method{logLik}{Mclust}(object, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'Mclust'} resulting from a call to \code{\link{Mclust}}.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{Returns an object of class \code{'logLik'} with an element providing the maximized log-likelihood, and further arguments giving the number of (estimated) parameters in the model (\code{"df"}) and the sample size (\code{"nobs"}).}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{Mclust}}.}

\examples{
\dontrun{
irisMclust <- Mclust(iris[,1:4])
summary(irisMclust)
logLik(irisMclust)
}
}

\keyword{multivariate}
