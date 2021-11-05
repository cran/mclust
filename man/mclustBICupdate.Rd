\name{mclustBICupdate}
\alias{mclustBICupdate}

\title{Update BIC values for parameterized Gaussian mixture models}

\description{
Update the BIC (Bayesian Information Criterion) for parameterized Gaussian 
mixture models by taking the best from BIC results as returned by \code{\link{mclustBIC}}.
}

\usage{
mclustBICupdate(BIC, \dots)
}

\arguments{
  \item{BIC}{Object of class \code{'mclustBIC'} containing the 
  BIC values as returned by a call to \code{\link{mclustBIC}}.
  }
  \item{\dots}{Further objects of class \code{'mclustBIC'} to be merged.}
}

\value{
An object of class \code{'mclustBIC'} containing the best values obtained from
merging the input arguments. Attributes are also updated according to the best
BIC found, so calling \code{\link{Mclust}} on the resulting ouput will return
the corresponding best model (see example).
}

\seealso{
  \code{\link{mclustBIC}},
  \code{\link{Mclust}}.
}
\examples{
\donttest{
data(galaxies, package = "MASS") 
galaxies <- galaxies / 1000

# use several random starting points
BIC <- NULL
for(j in 1:100)
{
  rBIC <- mclustBIC(galaxies, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BIC <- mclustBICupdate(BIC, rBIC)
}
pickBIC(BIC)
plot(BIC)

mod <- Mclust(galaxies, x = BIC)
summary(mod)
}
}
\keyword{cluster}
