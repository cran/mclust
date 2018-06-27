\name{plot.densityMclust}
\alias{plot.densityMclust}
\alias{plotDensityMclust1}
\alias{plotDensityMclust2}
\alias{plotDensityMclustd}

\title{Plots for Mixture-Based Density Estimate}

\description{
  Plotting methods for an object of class \code{'mclustDensity'}. Available graphs 
  are plot of BIC values and density for univariate and bivariate data. For 
  higher data dimensionality a scatterplot matrix of pairwise densities is
  drawn.
}

\usage{
\method{plot}{densityMclust}(x, data = NULL, what = c("BIC", "density", "diagnostic"), \dots)

plotDensityMclust1(x, data = NULL, hist.col = "lightgrey", 
                   hist.border = "white",  breaks = "Sturges", \dots)

plotDensityMclust2(x, data = NULL, nlevels = 11, levels = NULL, 
                   prob = c(0.25, 0.5, 0.75),
                   points.pch = 1, points.col = 1, points.cex = 0.8, \dots)

plotDensityMclustd(x, data = NULL, nlevels = 11, levels = NULL, 
                   prob = c(0.25, 0.5, 0.75),
                   points.pch = 1, points.col = 1, points.cex = 0.8,
                   gap = 0.2, \dots)
}

\arguments{
  \item{x}{An object of class \code{'mclustDensity'} obtained from a call to
           \code{\link{densityMclust}} function.}

  \item{data}{Optional data points.}

  \item{what}{The type of graph requested:
    \describe{
    \item{\code{"density"} =}{a plot of estimated density; if \code{data} is 
          also provided the density is plotted over data points (see Details 
          section).}
    \item{\code{"BIC"} =}{a plot of BIC values for the estimated models versus
                          the number of components.}
    \item{\code{"diagnostic"} =}{diagnostic plots (only available for the 
          one-dimensional case, see \code{\link{densityMclust.diagnostic}})}
    }
  }
  
  \item{hist.col}{The color to be used to fill the bars of the histogram.}

  \item{hist.border}{The color of the border around the bars of the histogram.}
  
  \item{breaks}{See the argument in function \code{\link[graphics]{hist}}.}

  \item{points.pch, points.col, points.cex}{The character symbols, colors, and magnification to be used for plotting \code{data} points.}

  \item{nlevels}{An integer, the number of levels to be used in plotting contour densities.}
  
  \item{levels}{A vector of density levels at which to draw the contour lines.}
  
  \item{prob}{A vector of probability levels for computing HDR. Only used if \code{type = "level"} and supersede previous \code{nlevels} and \code{levels} arguments.}

  \item{gap}{Distance between subplots, in margin lines, for the matrix of pairwise scatterplots.}

  \item{\dots}{Additional arguments passed to \code{\link{surfacePlot}}.}
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
  represented through \code{"contour"}, \code{"level"}, \code{"image"}, and 
  \code{"persp"} type of graph. 
  For \code{type = "level"} Highest Density Regions (HDRs) are plotted for 
  probability levels \code{prob}. See \code{link{hdrlevels}} for details.
  
  For higher dimensionality a scatterplot matrix of pairwise densities is
  drawn. 
}

% \value{}

\author{Luca Scrucca}

\seealso{
  \code{\link{densityMclust}},
  \code{\link{surfacePlot}},
  \code{\link{densityMclust.diagnostic}}, 
  \code{\link{Mclust}}.
}

\examples{
\dontrun{
dens <- densityMclust(faithful$waiting)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC", legendArgs = list(x = "topright"))
plot(dens, what = "density", data = faithful$waiting)

dens <- densityMclust(faithful)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "density", data = faithful, 
     drawlabels = FALSE, points.pch = 20)
plot(dens, what = "density", type = "level")
plot(dens, what = "density", type = "level", prob = seq(0.1, 0.9, by = 0.1))
plot(dens, what = "density", type = "level", data = faithful)
plot(dens, what = "density", type = "persp")

dens <- densityMclust(iris[,1:4])
summary(dens, parameters = TRUE)
plot(dens, what = "density", data = iris[,1:4], 
     col = "slategrey", drawlabels = FALSE, nlevels = 7)
plot(dens, what = "density", type = "level", data = iris[,1:4])
plot(dens, what = "density", type = "persp", col = grey(0.9))
}
}

\keyword{cluster}
\keyword{dplot}
