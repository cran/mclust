\name{plot.Mclust}
\alias{plot.Mclust}

\title{Plotting method for Mclust model-based clustering}

\description{
Plots for model-based clustering results, such as BIC, classification, uncertainty and density.
}
\usage{
\method{plot}{Mclust}(x, what = c("BIC", "classification", "uncertainty", "density"), 
     dimens = NULL, xlab = NULL, ylab = NULL,
     addEllipses = TRUE, main = FALSE, \dots)
}
\arguments{
  \item{x}{
    Output from \code{Mclust}.
  }
 \item{what}{
    A string specifying the type of graph requested. Available choices are:
    \describe{
      \item{\code{"BIC"}}{plot of BIC values used for choosing the number of clusters.}
      \item{\code{"classification"} =}{a plot showing the clustering. For data in more than two dimensions a pairs plot is produced, followed by a coordinate projection plot using specified \code{dimens}. Ellipses corresponding to covariances of mixture components are also drawn if \code{addEllipses = TRUE}.}
      \item{\code{"uncertainty"}}{a plot of classification uncertainty. For data in more than two dimensions a coordinate projection plot is drawn using specified \code{dimens}.}
      \item{\code{"density"}}{a plot of estimated density. For data in more than two dimensions a matrix of contours for coordinate projection plot is drawn using specified \code{dimens}.}
    }
    If not specified, in interactive sessions a menu of choices is proposed.
  }
\item{dimens}{
    A vector of integers specifying the dimensions of the coordinate projections
    in case of \code{"classification"}, \code{"uncertainty"}, or \code{"density"}
    plots. 
  }
 \item{xlab, ylab}{
    Optional labels for the x-axis and the y-axis.
  }
 \item{addEllipses}{
    A logical indicating whether or not to add ellipses with axes 
    corresponding to the within-cluster covariances in case of 
    \code{"classification"} or \code{"uncertainty"} plots. 
  }
 \item{main}{
    A logical or \code{NULL} indicating whether or not to add a title 
    to the plot identifying the type of plot drawn.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\details{
 For more flexibility in plotting, use \code{mclust1Dplot}, 
 \code{mclust2Dplot}, \code{surfacePlot}, \code{coordProj}, or
 \code{randProj}. 
}
\seealso{
  \code{\link{Mclust}},
  \code{\link{plot.mclustBIC}},
  \code{\link{plot.mclustICL}},
  \code{\link{mclust1Dplot}},
  \code{\link{mclust2Dplot}},
  \code{\link{surfacePlot}},
  \code{\link{coordProj}},
  \code{\link{randProj}}.
}
\examples{
\donttest{
precipMclust <- Mclust(precip)
plot(precipMclust)

faithfulMclust <- Mclust(faithful)
plot(faithfulMclust)

irisMclust <- Mclust(iris[,-5])
plot(irisMclust)
}
}
\keyword{cluster}
