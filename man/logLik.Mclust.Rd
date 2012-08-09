\name{logLik.Mclust}
\alias{logLik.Mclust}

\title{Log-Likelihood of a \code{Mclust} object}

\description{
Returns the log-likelihood for a \code{Mclust} object.}

\usage{
  \method{logLik}{Mclust}(object, data, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{"Mclust"} resulting from a call to \code{\link{Mclust}}.}

  \item{data}{the data for which the log-likelihood must be computed. If missing, the observed data from the \code{"Mclust"} object is used.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{Returns an object of class \code{logLik} with an element providing the maiximized log-likelihood, and further arguments giving the number of (estimated) parameters in the model (\code{"df"}) and the sample size (\code{"nobs"}).}

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

\seealso{\code{\link{Mclust}}.}

\examples{
\dontrun{
irisMclust <- Mclust(iris[,1:4])
summary(irisMclust)
logLik(irisMclust)
}
}

\keyword{multivariate}
