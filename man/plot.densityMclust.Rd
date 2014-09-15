\name{plot.densityMclust}
\alias{plot.densityMclust}
\alias{plotDensityMclust1}
\alias{plotDensityMclust2}
\alias{plotDensityMclustd}

\title{Plot for a \code{mclustDensity} object}

\description{
  Plotting methods for \code{mclustDensity} objects. Available graphs are
  plot of BIC values and density for univariate and bivariate data. For 
  higher data dimensionality a scatterplot matrix of pairwise densities is
  drawn.
}

\usage{
\method{plot}{densityMclust}(x, data = NULL, what = c("density", "BIC", "diagnostic"), \dots)

plotDensityMclust1(x, data = NULL, hist.col = "lightgrey", 
                   hist.border = "white",  breaks = "Sturges", \dots)

plotDensityMclust2(x, data = NULL, col = grey(0.6), nlevels = 11, 
                   levels = NULL, points.col = 1, pch = 1, \dots)

plotDensityMclustd(x, data = NULL, col = grey(0.6), nlevels = 11, 
                  levels = NULL, points.col = 1, pch = 1, gap = 0.2, \dots) 

}

\arguments{
  \item{x}{\code{mclustDensity} object obtained from \code{\link{densityMclust}} function.}

  \item{data}{optional data points.}

  \item{what}{the type of graph requested:
    \describe{
    \item{\code{"density"} =}{a plot of estimated density; if \code{data} is also provided the density is plotted over data points (see Details section).}
    \item{\code{"BIC"} =}{a plot of BIC values for the estimated models versus the number of components.}
    \item{\code{"diagnostic"} =}{diagnostic plots (only available for the one-dimensional case, see \code{\link{densityMclust.diagnostic}})}
    }
  }
  
  \item{hist.col}{the colour to be used to fill the bars of the histogram.}

  \item{hist.border}{the color of the border around the bars of the histogram.}
  
  \item{breaks}{see the argument in function \code{\link[graphics]{hist}}.}

  \item{points.col}{the color to be used for plotting \code{data} points.}

  \item{pch}{an integer, a symbol or a single character to be used for plotting \code{data} points.}

  \item{col}{color for the contour lines drawn.} 
  
  \item{nlevels}{an integer, the number of levels to be used in plotting contour densities.}
  
  \item{levels}{a vector of density levels at which to draw the contour lines.}
  
  \item{gap}{distance between subplots, in margin lines, for the matrix of pairwise scatterplots.}

  \item{\dots}{additional arguments.}
}

\details{The function \code{plot.densityMclust} allows to obtain the plot of
  estimated density or the graph of BIC values for evaluated models. 
  
  If \code{what = "density"} the produced plot dependes on the dimensionality
  of the data.
  
  For one-dimensional data a call with no \code{data} provided produces a 
  plot of the estimated density over a sensible range of values. If 
  \code{data} is provided the density is over-plotted on a histogram for the
  observed data. 
  
  For two-dimensional data further arguments available are those accepted by
  the \code{\link{surfacePlot}} function. In particular, the density can be
  represented through \code{"contour"}, \code{"image"}, and \code{"persp"}
  type of graph.
  
  For higher dimensionality a scatterplot matrix of pairwise densities is
  drawn. 
}

% \value{}

\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Luca Scrucca}

\seealso{
  \code{\link{densityMclust}},
  \code{\link{densityMclust.diagnostic}}, 
  \code{\link{Mclust}}.
}

\examples{
x = faithful$waiting
dens = densityMclust(x)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens)
plot(dens, x)
plot(dens, what = "diagnostic")

x = as.matrix(faithful)
dens = densityMclust(x)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens)
plot(dens, x, col = "cadetblue", drawlabels = FALSE, pch = 20,
     levels = quantile(dens$density, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
plot(dens, x, col = "grey", 
     points.col = dens$classification+1, pch = dens$classification)
plot(dens, type = "image", col = terrain.colors(50))
plot(dens, type = "persp")

x = iris[,1:4]
dens = densityMclust(x)
summary(dens, parameters = TRUE)
plot(dens)
plot(dens, x, col = "cadetblue", drawlabels = FALSE,
     levels = quantile(dens$density, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
plot(dens, type = "image", col = topo.colors(50))
plot(dens, type = "persp", border = NA)
}

\keyword{cluster}
\keyword{dplot}
