\name{mclust2Dplot}
\alias{mclust2Dplot}

\title{Plot two-dimensional data modelled by an MVN mixture}

\description{
  Plot two-dimensional data given parameters of an MVN mixture model 
  for the data.
}
\usage{
mclust2Dplot(data, parameters = NULL, z = NULL,
             classification = NULL, truth = NULL, uncertainty = NULL,
             what = c("classification", "uncertainty", "error"), 
             addEllipses = TRUE, fillEllipses = mclust.options("fillEllipses"),
             symbols = NULL, colors = NULL, 
             xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
             scale = FALSE, cex  = 1, PCH = ".",
             main = FALSE, swapAxes = FALSE, \dots) 
}
\arguments{
  \item{data}{
    A numeric matrix or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables. 
    In this case the data are two dimensional, so there are two columns.
  }
  \item{parameters}{
     A named list giving the parameters of an \emph{MCLUST} model, 
     used to produce superimposing ellipses on the plot. 
     The relevant components are as follows:
     \describe{
        \item{\code{pro}}{
              Mixing proportions for the components of the mixture. 
              There should one more mixing proportion than the number of 
              Gaussian components if the mixture model includes 
              a Poisson noise term.
        }
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
 \item{z}{
        A matrix in which the \code{[i,k]}th entry gives the
        probability of observation \emph{i} belonging to the \emph{k}th class. 
        Used to compute \code{classification} and
        \code{uncertainty} if those arguments aren't available.
  }
 \item{classification}{
        A numeric or character vector representing a classification of
        observations (rows) of \code{data}. If present argument \code{z}
        will be ignored.
  }
\item{truth}{
        A numeric or character vector giving a known
        classification of each data point.
        If \code{classification}
        or \code{z} is also present, 
        this is used for displaying classification errors.
  }
  \item{uncertainty}{
        A numeric vector of values in \emph{(0,1)} giving the
        uncertainty of each data point. If present argument \code{z}
        will be ignored.
  }
  \item{what}{
    Choose from one of the following three options: \code{"classification"}
    (default), \code{"error"}, \code{"uncertainty"}. 
  }
  \item{addEllipses}{
    A logical indicating whether or not to add ellipses with axes 
    corresponding to the within-cluster covariances.
  }
  \item{fillEllipses}{
    A logical specifying whether or not to fill ellipses with transparent
    colors when \code{addEllipses = TRUE}.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class in \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{unique}). 
    The default is given by \code{mclust.options("classPlotSymbols")}.
  }
  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique class in \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{unique}). 
    The default is given is \code{mclust.options("classPlotColors")}.
  }
  \item{xlim, ylim}{
    Optional argument specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{xlab, ylab}{
    Optional argument specifying labels for the x-axis and y-axis.
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=FALSE} 
  }
  \item{cex}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{PCH}{
    An argument specifying the symbol to be used when a classificatiion
    has not been specified for the data. The default value is a small dot ".".
  }
  \item{main}{
    A logical variable or \code{NULL} indicating whether or not to add a title 
    to the plot identifying the dimensions used.
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
A plot showing the data, together with the location of the mixture components, classification, uncertainty, and/or classification errors.
}
\seealso{
  \code{\link{surfacePlot}},
  \code{\link{clPairs}},
  \code{\link{coordProj}},
  \code{\link{mclust.options}}
}
\examples{
\donttest{
faithfulModel <- Mclust(faithful)

mclust2Dplot(faithful, parameters=faithfulModel$parameters, 
             z=faithfulModel$z, what = "classification", main = TRUE)

mclust2Dplot(faithful, parameters=faithfulModel$parameters, 
             z=faithfulModel$z, what = "uncertainty", main = TRUE)
}
}
\keyword{cluster}
