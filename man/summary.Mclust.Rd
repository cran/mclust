\name{summary.Mclust}
\alias{summary.Mclust}
\alias{print.summary.Mclust}

\title{Summarizing Gaussian Finite Mixture Model Fits}

\description{Summary method for class \code{"Mclust"}.}

\usage{
\method{summary}{Mclust}(object, parameters = FALSE, classification = FALSE, \dots)
\method{print}{summary.Mclust}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{"Mclust"} resulting of a call to \code{\link{Mclust}} or \code{\link{densityMclust}}.}

  \item{x}{An object of class \code{"summary.Mclust"}, usually, a result of a call to \code{summary.Mclust}.}

  \item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

  \item{classification}{Logical; if \code{TRUE}, the MAP classification/clustering of observations is printed.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

% \details{}

% \value{}

\references{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{Mclust}}, \code{\link{densityMclust}}.}

\examples{
mod1 = Mclust(iris[,1:4])
summary(mod1)
summary(mod1, parameters = TRUE, classification = TRUE)

mod2 = Mclust(iris[,1:4], G = 1)
summary(mod2, parameters = TRUE, classification = TRUE)

mod3 = Mclust(iris[,1:4], prior = priorControl())
summary(mod3)

mod4 = Mclust(iris[,1:4], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(mod4, parameters = TRUE, classification = TRUE)
}

\keyword{cluster}
