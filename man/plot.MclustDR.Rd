\name{plot.MclustDR}
\alias{plot.MclustDR}
\alias{plotEvalues.MclustDR}

\title{Plotting method for dimension reduction for model-based clustering and classification}

\description{
Graphs data projected onto the estimated subspace for model-based clustering and classification.
}
\usage{
\method{plot}{MclustDR}(x, dimens, 
     what = c("scatterplot", "pairs", "contour", "classification",
              "boundaries", "density", "evalues"), 
     symbols, colors, col.contour = gray(0.7), col.sep = grey(0.4), 
     ngrid = 100, nlevels = 5, asp = NULL, \dots)
}

\arguments{
  \item{x}{
    An object of class \code{'MclustDR'} resulting from a call to \code{\link{MclustDR}}.
}
  \item{dimens}{
    A vector of integers giving the dimensions of the desired coordinate
    projections for multivariate data.
}
  \item{what}{
        The type of graph requested: 
    \describe{
      \item{\code{"scatterplot"} =}{a two-dimensional plot of data projected onto the first two directions specified by \code{dimens} and with data points marked according to the corresponding mixture component.
      By default, the first two directions are selected for plotting.}
      \item{\code{"pairs"} =}{a scatterplot matrix of data projected onto the estimated subspace and with data points marked according to the corresponding mixture component.
         By default, all the available directions are used, unless they have been specified by \code{dimens}.}
      \item{\code{"contour"} =}{a two-dimensional plot of data projected onto the first two directions specified by \code{dimens} (by default, the first two directions) with density contours for classes or clusters and data points marked according to the corresponding mixture component.}
      \item{\code{"classification"} =}{a two-dimensional plot of data projected onto the first two directions specified by \code{dimens} (by default, the first two directions) with classification region and data points marked according to the corresponding mixture component.}
      \item{\code{"boundaries"} =}{a two-dimensional plot of data projected onto the first two directions specified by \code{dimens} (by default, the first two directions) with uncertainty boundaries and data points marked according to the corresponding mixture component.
      The uncertainty is shown using a greyscale with darker regions indicating higher uncertainty.
      }
      \item{\code{"density"} =}{a one-dimensional plot of estimated density for the first direction specified by \code{dimens} (by default, the first one). A set of box-plots for each estimated cluster or known class are also shown at the bottom of the graph.
      }
}
}
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique mixture component. Elements in \code{colors} correspond to classes
    in order of appearance in the sequence of observations (the order used by
    the function \code{factor}). 
    The default is given by \code{mclust.options("classPlotSymbols")}.
}
  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique cluster or known class. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{factor}). 
    The default is given by \code{mclust.options("classPlotColors")}.
}
  \item{col.contour}{
    The color of contours in case \code{what = "contour"}.
}
  \item{col.sep}{
    The color of classification boundaries in case \code{what = "classification"}.
}
  \item{ngrid}{
    An integer specifying the number of grid points to use in evaluating the classification regions.
}
  \item{nlevels}{
    The number of levels to use in case \code{what = "contour"}.
}
  \item{asp}{For scatterplots the \eqn{y/x} aspect ratio, see
    \code{\link{plot.window}}.
}
  \item{\dots}{further arguments passed to or from other methods.}
}

%\details{}

%\value{}

\references{
  Scrucca, L. (2010) Dimension reduction for model-based clustering.
  \emph{Statistics and Computing}, 20(4), pp. 471-484.

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Luca Scrucca}

%\note{}

\seealso{\link{MclustDR}}

\examples{
mod = Mclust(iris[,1:4], G = 1:5)
dr = MclustDR(mod)
plot(dr, what = "evalues")
plot(dr, what = "scatterplot", dimens = c(1,3))
plot(dr, what = "pairs")
plot(dr, what = "contour")
plot(dr, what = "classification", ngrid = 200)
plot(dr, what = "boundaries", ngrid = 200)
plot(dr, what = "density")
plot(dr, what = "density", dimens = 2)

data(banknote)
da = MclustDA(banknote[,2:7], banknote$Status, G = 1:3)
dr = MclustDR(da)
plot(dr, what = "evalues")
plot(dr, what = "pairs")
plot(dr, what = "contour")
plot(dr, what = "contour", dimens = c(1,3))
plot(dr, what = "classification", ngrid = 200)
plot(dr, what = "boundaries", ngrid = 200)
plot(dr, what = "density")
plot(dr, what = "density", dimens = 2)
}

\keyword{multivariate}
