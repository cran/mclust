\name{densityMclust}
\alias{densityMclust}

\title{Density Estimation via Model-Based Clustering}

\description{
  Produces a density estimate for each data point using a Gaussian finite 
  mixture model from \code{Mclust}.
}

\usage{
densityMclust(data, \dots)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{\dots }{
    Additional arguments for the \code{\link{Mclust}} function. 
    In particular, setting the arguments \code{G} and \code{modelNames} 
    allow to specify the number of mixture components and the type of
    model to be fitted. By default an "optimal" model is selected based
    on the BIC criterion. 
  }
}

\value{
  An object of class \code{densityMclust}, which inherits from 
  \code{Mclust}, is returned with the following slots added:
  \item{varname}{
    A character string denoting the variables names.}
  \item{range}{
     The range of the input data (used as a default for plotting).}
  \item{density}{
    The density of the input data \code{x} according to the estimated 
    model.}
}

%\details{}
  
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Revised version by Luca Scrucca based on 
  the original code by C. Fraley and A.E. Raftery.}

\seealso{
  \code{\link{plot.densityMclust}}, 
  \code{\link{Mclust}}, 
  \code{\link{summary.Mclust}},
  \code{\link{predict.densityMclust}}.
}

\examples{
x = faithful$waiting
dens = densityMclust(x)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens)
plot(dens, x)

x = as.matrix(faithful)
dens = densityMclust(x)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens)
plot(dens, x, col = "cadetblue", drawlabels = FALSE, pch = 20,
     levels = quantile(dens$density, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
plot(dens, x, col = "grey", 
     points.col = mclust.options("classPlotColors")[dens$classification], 
     pch = dens$classification)
plot(dens, type = "image", col = topo.colors(50))
plot(dens, type = "persp")

x = iris[,1:4]
dens = densityMclust(x)
summary(dens, parameters = TRUE)
plot(dens)
plot(dens, x, col = "cadetblue", drawlabels = FALSE,
     levels = quantile(dens$density, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
plot(dens, type = "image", col = gray.colors(50))
plot(dens, type = "persp", col = gray(0.5))
}

\keyword{cluster}
