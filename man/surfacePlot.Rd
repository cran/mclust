\name{surfacePlot}
\alias{surfacePlot}

\title{Density or uncertainty surface for bivariate mixtures}

\description{
Plots a density or uncertainty surface given bivariate data and parameters of 
a MVN mixture model for the data.
}
\usage{
surfacePlot(data, parameters, 
            type = c("contour", "hdr", "image", "persp"), 
            what = c("density", "uncertainty"), 
            transformation = c("none", "log", "sqrt"),          
            grid = 100, nlevels = 11, levels = NULL, 
            prob = c(0.25, 0.5, 0.75),
            color.palette = blue2grey.colors, col = grey(0.6),
            xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
            scale = FALSE, main = FALSE, swapAxes = FALSE, 
            verbose = FALSE, \dots)
}
\arguments{
  \item{data}{
    A matrix, or data frame of bivariate observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
 \item{parameters}{
     A named list giving the parameters of an \emph{MCLUST} model, 
     used to produce superimposing ellipses on the plot. 
     The relevant components are as follows:
     \describe{
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}}
              for details.
        }
     }
  }
  \item{type}{
    Choose from one of the following three options: \code{"contour"} 
    (default), \code{"hdr"}, \code{"image"}, and \code{"persp"} indicating 
    the plot type. 
  }
  \item{what}{
    Choose from one of the following options: \code{"density"}
    (default), \code{"uncertainty"} indicating what to plot. 
  }
  \item{transformation}{
    Choose from one of the following three options: \code{"none"}
    (default), \code{"log"}, \code{"sqrt"} indicating a transformation
    to be applied before plotting. 
  }
  \item{grid}{
    The number of grid points (evenly spaced on each axis). 
    The mixture density and uncertainty is computed at 
    \code{grid x grid} points to produce the surface plot.
    Default: \code{100}.
  }
  \item{nlevels}{
    The number of levels to use for a contour plot.
    Default: \code{11}.
  }
  \item{levels}{
    A vector of levels at which to draw the lines in a contour plot.
  }
  \item{prob}{
    A vector of probability levels for computing HDR. 
    Only used if \code{type = "hdr"} and supersede previous 
    \code{nlevels} and \code{levels} arguments.
  }
  \item{color.palette}{
    A function which defines a palette of colours to be used for 
    \code{type = "hdr"} plot. 
  }
  \item{col}{
    The color to be used for drawing contour lines.
  }
  \item{xlim, ylim}{
    Optional argument specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{xlab, ylab}{
    Optional argument specifying labels for the x-axis and y-axis.
  }
  \item{scale}{
    A logical variable indicating whether or not the two
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    The default is not to scale.
  }
  \item{main}{
    A logical variable or \code{NULL} indicating whether or not to add a title 
    to the plot identifying the dimensions used.
  }
  \item{swapAxes}{
    A logical variable indicating whether or not the axes should be swapped
    for the plot.
  }
  \item{verbose}{
    A logical variable telling whether or not to print an indication that
    the function is in the process of computing values at the grid points,
    which typically takes some time to complete.
  }
 \item{\dots}{
    Other graphics parameters.
  }
}
\value{
A plots showing (a transformation of) the density or uncertainty for the given
mixture model and data. 

The function also returns an invisible list with components \code{x}, 
\code{y}, and \code{z} in which \code{x} and \code{y} are the values used to 
define the grid and \code{z} is the transformed density or uncertainty at the 
grid points.
}
\details{
  For an image plot, a color scheme may need to be selected on the display
  device in order to view the plot.
}
\references{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\seealso{
  \code{\link{mclust2Dplot}}
}
\examples{
\dontrun{
faithfulModel <- Mclust(faithful)
surfacePlot(faithful, parameters = faithfulModel$parameters,
            type = "contour", what = "density", transformation = "none",
            drawlabels = FALSE)
surfacePlot(faithful, parameters = faithfulModel$parameters,
            type = "persp", what = "density", transformation = "log")
surfacePlot(faithful, parameters = faithfulModel$parameters,
            type = "contour", what = "uncertainty", transformation = "log")
}
}
\keyword{cluster}
