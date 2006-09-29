\name{surfacePlot}
\alias{surfacePlot}
\title{
  Density or uncertainty surface for two dimensional mixtures.
}
\description{
  Plots a density or uncertainty surface 
  given data in more than two dimensions
  and parameters of an MVN mixture model for the data.
}
\usage{
surfacePlot(data, parameters, 
            type = c("contour", "image", "persp"), 
            what = c("density", "uncertainty"), 
            transformation = c("none", "log", "sqrt"),          
            grid = 50, nlevels = 20, scale = FALSE, 
            xlim=NULL, ylim=NULL, 
            identify = FALSE, verbose = FALSE, swapAxes = FALSE, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
 \item{parameters}{
     A named list giving the parameters of an \emph{MCLUST} model, 
     used to produce superimposing ellipses on the plot. 
     The relevant components are as follows:
     \describe{
        \item{mean}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{variance}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{mclustVariance}
              for details.
        }
     }
  }
  \item{type}{
    Choose from one of the following three options: \code{"contour"}
    (default), \code{"image"}, \code{"persp"} indicating the plot type. 
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
    Default: \code{50}.
  }
  \item{nlevels}{
    The number of levels to use for a contour plot.
    Default: \code{20}.
  }
  \item{scale}{
    A logical variable indicating whether or not the two
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    The default is not to scale.
  }
  \item{xlim, ylim}{
    An argument specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{identify}{
    A logical variable indicating whether or not to add a title to the plot
    identifying the dimensions used.
  }
  \item{verbose}{
    A logical variable telling whether or not to print an indication that
    the function is in the process of computing values at the grid points,
    which typically takes some time to complete.
  }
  \item{swapAxes}{
    A logical variable indicating whether or not the axes should be swapped
    for the plot.
  }
 \item{\dots}{
    Other graphics parameters.
  }
}
\value{
  An invisible list with components x, y, and z in which x and y are the
  values used to define the grid and z is the transformed density or
  uncertainty at the grid points.
}
\section{Side Effects}{
  A plots showing (a transformation of) the density or uncertainty
  for the given mixture model and data. 
}
\section{Details}{
  For an image plot, a color scheme may need to be selected on the display
  device in order to view the plot.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{mclust2Dplot}}
}
\examples{
faithfulModel <- mclustModel(faithful,mclustBIC(faithful))
surfacePlot(faithful, parameters = faithfulModel$parameters,
            type = "contour", what = "density", transformation = "none",
            drawlabels = FALSE)

}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
