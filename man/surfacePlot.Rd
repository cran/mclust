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
surfacePlot(data, mu, pro, \dots, type = c("contour", "image", "persp"), 
        what = c("density", "uncertainty", "skip"), 
        transformation = c("none", "log", "sqrt"),          
        grid = 50, nlevels = 20, scale = FALSE, identify = FALSE,
        verbose = TRUE, xlim, ylim, swapAxes = FALSE)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{mu}{
    A matrix whose columns are the means of each group. 
  }
  \item{pro}{
    A list with \code{scale}, \code{shape} and \code{orientation}
    components giving an alternative form for the covariance structure
    of the mixture model.
  }
  \item{\dots}{
    An argument specifying the covariance structure of the model. If
    used an indirect function call via \code{do.call} (see example
    below), it is usually not necessary to know the precise form for
    this argument. This argument usually take one of the following
    forms:
    \describe{
      \item{\code{sigma}}{ A three dimensional array in which
	\code{sigma[,,k]} gives the covariance for the \emph{k}th group.}
      \item{\code{decomp}}{ A list with \code{scale}, \code{shape} and
	\code{orientation} components giving an alternative form for the
	covariance structure of the mixture model.}
    }
  }
  \item{type}{
    Any subset of \code{c("contour","image","persp")} indicating the
    plot type. For more than one selection, users will be asked to
    choose from a menu.
  }
  \item{what}{
    Any subset of \code{c("density","uncertainty","skip")} indicating
    what to plot. For more than one selection, users will be asked to
    choose from a menu. The \code{"skip"} produces and empty plot, which
    may be useful if multiple plots are displayed simultaneously. 
  }
  \item{transformation}{
    Any subset of \code{c("none","log","sqrt")} indicating  a
    transformation to be applied to the surface values before
    plotting. For more than one selection, users will be asked to choose
    from a menu. 
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
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=F} 
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
  \item{xlim, ylim }{
    An argument specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{swapAxes}{
    A logical variable indicating whether or not the axes should be swapped
    for the plot.
  }
}
\value{
  An invisible list with components x, y, and z in which x and y are the
  values used to define the grid and z is the transformed density or
  uncertainty at the grid points.
}
\section{Side Effects}{
  One or more plots showing location of the mixture components, classification,
  uncertainty, and/or classification errors.
}
\section{Details}{
  For an image plot, a color scheme may need to be selected on the display
  device in order to view the plot.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}.
  
  
  C. Fraley and A. E. Raftery (2002).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis.
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{mclust2Dplot}},
  \code{\link{do.call}}
}
\examples{
n <- 250 ## create artificial data
set.seed(0)
x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))

xEMclust <- summary(EMclust(x),x)
surfacePlot(x, mu = xEMclust$mu, sigma = xEMclust$sigma, pro=xEMclust$pro,
            type = "contour", what = "density", transformation = "none")

\dontrun{do.call("surfacePlot", c(list(data = x), xEMclust))}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
