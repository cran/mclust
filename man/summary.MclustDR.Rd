\name{summary.MclustDR}
\alias{summary.MclustDR}
\alias{print.summary.MclustDR}

\title{Summarizing dimension reduction method for model-based clustering and classification}

\description{Summary method for class \code{"MclustDR"}.}

\usage{
\method{summary}{MclustDR}(object, numdir, std = FALSE, \dots)
\method{print}{summary.MclustDR}(x, digits = max(5, getOption("digits") - 3), \dots)
}

\arguments{
  \item{object}{An object of class \code{'MclustDR'} resulting from a call to \code{\link{MclustDR}}.}

  \item{x}{An object of class \code{'summary.MclustDR'}, usually, a result of a call to \code{summary.MclustDR}.}

  \item{numdir}{An integer providing the number of basis directions to be printed.}

  \item{std}{if \code{TRUE} the coefficients basis are scaled such that all
             predictors have unit standard deviation.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

%\details{}
%\value{}
\author{Luca Scrucca}
%\note{}

\seealso{
\code{\link{MclustDR}}, \code{\link{plot.MclustDR}}
}
%\examples{}
%\keyword{}
