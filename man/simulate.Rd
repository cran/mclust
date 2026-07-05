\name{simulate.Mclust}
\alias{simulate.Mclust}
\alias{simulate.densityMclust}
\title{
Simulate from Gaussian mixture models
}
\description{
Generic method for simulating from Gaussian mixture models via base `simulate()` function.
}
\usage{
\method{simulate}{Mclust}(object, nsim = 1, seed = NULL, \dots)
\method{simulate}{densityMclust}(object, nsim = 1, seed = NULL, \dots)
}
\arguments{
  \item{object}{An object of either `Mclust` or `densityMclust` class.}
 \item{nsim}{An integer specifying the number of data points to be simulated.}
  \item{seed}{An optional integer argument to \code{set.seed} for reproducible
    random class assignment.  By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{sim}.
  }
  \item{\dots }{Catches unused arguments in indirect or list calls via \code{do.call}.}
}
\value{
  For an input object of class `Mclust` returns a matrix in which the 
  first column is the mixture component label and the remaining
  columns are the simulated data observations. 
  For an input object of class `densityMclust` only the simulated data 
  observations are returned.
}

\seealso{
  \code{\link{simulate}},
  \code{\link{sim}}, \dots,
  \code{\link{Mclust}}, 
  \code{\link{densityMclust}}
}

\examples{
mod <- Mclust(iris[,-5])
iris_sim <- simulate(mod, n = 1000)
str(iris_sim)
clPairs(iris_sim[,-1], iris_sim[,1])

x <- iris$Petal.Length
mod <- densityMclust(x, plot = FALSE)
x_sim <- simulate(mod, n = 1000)
str(x_sim)
hist(x_sim, breaks = 31, main = NULL, probability = TRUE)
rug(x_sim)
x0 = seq(0.5, 9, length = 101)
lines(x0, predict(mod, newdata = x0))
}
