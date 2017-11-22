\name{summary.MclustDA}
\alias{summary.MclustDA}
\alias{print.summary.MclustDA}

\title{Summarizing discriminant analysis based on Gaussian finite mixture modeling.}

\description{Summary method for class \code{"MclustDA"}.}

\usage{
\method{summary}{MclustDA}(object, parameters = FALSE, newdata, newclass, \dots)
\method{print}{summary.MclustDA}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.}

  \item{x}{An object of class \code{'summary.MclustDA'}, usually, a result of a call to \code{summary.MclustDA}.}

  \item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

  \item{newdata}{A data frame or matrix giving the test data.}

  \item{newclass}{A vector giving the class labels for the observations in 
    the test data.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

% \details{}

\value{The function \code{summary.MclustDA} computes and returns a list of summary statistics of the estimated MclustDA or EDDA model for classification.}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{MclustDA}}, \code{\link{plot.MclustDA}}.}

\examples{
mod = MclustDA(data = iris[,1:4], class = iris$Species)
summary(mod)
summary(mod, parameters = TRUE)
}
\keyword{multivariate}
