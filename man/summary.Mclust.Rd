\name{summary.Mclust}
\alias{summary.Mclust}
\alias{print.summary.Mclust}

\title{Summarizing Gaussian Finite Mixture Model Fits}

\description{Summary method for class \code{"Mclust"}.}

\usage{
\method{summary}{Mclust}(object, classification = TRUE, parameters = FALSE, \dots)
\method{print}{summary.Mclust}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{'Mclust'} resulting of a call to \code{\link{Mclust}} or \code{\link{densityMclust}}.}

  \item{x}{An object of class \code{'summary.Mclust'}, usually, a result of a call to \code{summary.Mclust}.}

  \item{classification}{Logical; if \code{TRUE} a table of MAP classification/clustering of observations is printed.}

  \item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

% \details{}

% \value{}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{Mclust}}, \code{\link{densityMclust}}.}

\examples{
mod1 = Mclust(iris[,1:4])
summary(mod1)
summary(mod1, parameters = TRUE, classification = FALSE)

mod2 = densityMclust(faithful)
summary(mod2)
summary(mod2, parameters = TRUE)
}

\keyword{cluster}
