\name{summary.MclustSSC}
\alias{summary.MclustSSC}
\alias{print.summary.MclustSSC}

\title{Summarizing semi-supervised classification model based on Gaussian finite mixtures}

\description{Summary method for class \code{"MclustSSC"}.}

\usage{
\method{summary}{MclustSSC}(object, parameters = FALSE, \dots)
\method{print}{summary.MclustSSC}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{'MclustSSC'} resulting from a call to \code{\link{MclustSSC}}.}

  \item{x}{An object of class \code{'summary.MclustSSC'}, usually, a result of a call to \code{summary.MclustSSC}.}

  \item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

% \details{}

\value{The function \code{summary.MclustSSC} computes and returns a list of summary statistics of the estimated MclustSSC model for semi-supervised classification.}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{MclustSSC}}, \code{\link{plot.MclustSSC}}.}

\keyword{classification}
