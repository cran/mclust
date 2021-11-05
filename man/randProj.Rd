\name{randProj}
\alias{randProj}

\title{Random projections of multidimensional data modeled by an MVN mixture}
  
\description{
  Plots random projections given multidimensional data
  and parameters of an MVN mixture model for the data.
}
\usage{
randProj(data, seeds = NULL, parameters = NULL, z = NULL,
         classification = NULL, truth = NULL, uncertainty = NULL, 
         what = c("classification", "error", "uncertainty"),
         quantiles = c(0.75, 0.95), 
         addEllipses = TRUE, fillEllipses = mclust.options("fillEllipses"),
         symbols = NULL, colors = NULL, scale = FALSE, 
         xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
         cex = 1, PCH = ".", main = FALSE, \dots)
}
\arguments{
  \item{data}{
    A numeric matrix or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{seeds}{
    An integer value or a vector of integer values to be used as seed for 
    random number generation. If multiple values are provided, then each seed 
    should produce a different projection. 
    By default, a single seed is drawn randomnly, so each call of 
    \code{randProj()} produces different projections.
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
  \item{quantiles}{
    A vector of length 2 giving quantiles used in plotting
    uncertainty. The smallest symbols correspond to the smallest
    quantile (lowest uncertainty), medium-sized (open) symbols to points
    falling between the given quantiles, and large (filled) symbols to
    those in the largest quantile (highest uncertainty). The default is
    \emph{(0.75,0.95)}. 
  }
  \item{addEllipses}{
    A logical indicating whether or not to add ellipses with axes 
    corresponding to the within-cluster covariances in case of 
    \code{"classification"} or \code{"uncertainty"} plots. 
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
    The default is given by \code{mclust.options("classPlotColors")}.
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=FALSE} 
  }
  \item{xlim, ylim}{
    Optional arguments specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{xlab, ylab}{
    Optional arguments specifying the labels for, respectively, the horizontal 
    and vertical axis.
  } 
  \item{cex}{
    A numerical value specifying the size of the plotting symbols. 
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
  \item{\dots}{
    Other graphics parameters.
  }
}
\value{
A plot showing a random two-dimensional projection of the data, together with the location of the  mixture components, classification, uncertainty, and/or classification errors. 

The function also returns an invisible list with components \code{basis}, the randomnly generated basis of the projection subspace, \code{data}, a matrix of projected data, and \code{mu} and \code{sigma} the component parameters transformed to the projection subspace.
}

\seealso{
  \code{\link{clPairs}},
  \code{\link{coordProj}},
  \code{\link{mclust2Dplot}},
  \code{\link{mclust.options}}
}
\examples{
\donttest{
est <- meVVV(iris[,-5], unmap(iris[,5]))
par(pty = "s", mfrow = c(1,1))
randProj(iris[,-5], seeds=1:3, parameters = est$parameters, z = est$z,
          what = "classification", main = TRUE) 
randProj(iris[,-5], seeds=1:3, parameters = est$parameters, z = est$z,
          truth = iris[,5], what = "error", main = TRUE) 
randProj(iris[,-5], seeds=1:3, parameters = est$parameters, z = est$z,
          what = "uncertainty", main = TRUE) 
}
}
\keyword{cluster}
