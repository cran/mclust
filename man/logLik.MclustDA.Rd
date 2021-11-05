\name{logLik.MclustDA}
\alias{logLik.MclustDA}

\title{Log-Likelihood of a \code{MclustDA} object}

\description{
Returns the log-likelihood for a \code{MclustDA} object.}

\usage{
  \method{logLik}{MclustDA}(object, data, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.}

  \item{data}{the data for which the log-likelihood must be computed. If missing, the observed data from the \code{'MclustDA'} object is used.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{Returns an object of class \code{'logLik'} with an element providing the maximized log-likelihood, and further arguments giving the number of (estimated) parameters in the model (\code{"df"}) and the sample size (\code{"nobs"}).}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{MclustDA}}.}

\examples{
\donttest{
irisMclustDA <- MclustDA(iris[,1:4], iris$Species)
summary(irisMclustDA)
logLik(irisMclustDA)
}
}

\keyword{multivariate}
