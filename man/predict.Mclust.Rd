\name{predict.Mclust}
\alias{predict.Mclust}

\title{Cluster multivariate observations by Gaussian finite mixture modeling}

\description{Cluster prediction for multivariate observations based on Gaussian finite mixture models estimated by \code{\link{Mclust}}.}

\usage{
  \method{predict}{Mclust}(object, newdata, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{"Mclust"} resulting from a call to \code{\link{Mclust}}.}

  \item{newdata}{a data frame or matrix giving the data. If missing the clustering data obtained from the call to \code{\link{Mclust}} are classified.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
Returns a list of with the following components:
  \item{classification}{a factor of predicted cluster labels for \code{newdata}.}
  \item{z}{a matrix whose \emph{[i,k]}th entry is the probability that 
           observation \emph{i} in \code{newdata} belongs to the \emph{k}th cluster.}
}

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
model <- Mclust(faithful)

# predict cluster for the observed data
pred <- predict(model) 
str(pred)
pred$z              # equal to model$z
pred$classification # equal to 
plot(faithful, col = pred$classification, pch = pred$classification)

# predict cluster over a grid
grid <- apply(faithful, 2, function(x) seq(min(x), max(x), length = 50))
grid <- expand.grid(eruptions = grid[,1], waiting = grid[,2])
pred <- predict(model, grid)
plot(grid, col = mclust.options()$classPlotColors[pred$classification], pch = 15, cex = 0.5)
points(faithful, pch = model$classification)
}

\keyword{multivariate}
