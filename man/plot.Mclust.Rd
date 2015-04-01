\name{plot.Mclust}
\alias{plot.Mclust}
\title{
  Plot Model-Based Clustering Results
}
\description{
  Plot model-based clustering results: BIC, classification, uncertainty and 
  (for univariate and bivariate data) density.
}
\usage{
\method{plot}{Mclust}(x, what = c("BIC", "classification", "uncertainty", "density"), 
     dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
     addEllipses = TRUE, main = TRUE, \dots)
}
\arguments{
  \item{x}{
    Output from \code{Mclust}.
  }
 \item{what}{
    The type of graph requested:
    \describe{
      \item{\code{"BIC"}}{}
      \item{\code{"classification"}}{}
      \item{\code{"uncertainty"}}{}
      \item{\code{"density"}}{}
    }
    By default, all the above graphs are produced. See the description below.
  }
\item{dimens}{
    A vector of length one or two giving the integer dimensions of the
    desired coordinate projections for multivariate data in case of
    \code{"classification"} or \code{"uncertainty"} plots. 
  }
 \item{xlab, ylab}{
    Optional labels for the x-axis and the y-axis.
  }
 \item{ylim}{
    Optional limits for the vertical axis of the BIC plot.
  }
 \item{addEllipses}{
    A logical indicating whether or not to add ellipses with axes corresponding 
    to the within-cluster covariances if \code{what = "classification"}. 
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
  Model-based clustering plots: 
  \describe{
  \item{\code{"BIC"} =}{BIC values used for choosing the number of clusters.}
  \item{\code{"classification"} =}{a plot showing the clustering. For data
      in more than two dimensions a pairs plot is produced, followed by a
      coordinate projection plot using specified \code{dimens}.} 
  \item{\code{"uncertainty"} =}{a plot of classification uncertainty. For
      data in more than two dimensions a coordinate projection plot is 
      drawn using specified \code{dimens}.}
  \item{\code{"density"} =}{a plot of estimated density. For two
      dimensional a contour plot is drawn, while for data in more than 
      two dimensions a matrix of contours for pairs of variables is 
      produced.}
    }
}
\details{
 For more flexibility in plotting, use \code{mclust1Dplot}, 
 \code{mclust2Dplot}, \code{surfacePlot}, \code{coordProj}, or
 \code{randProj}. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
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
precipMclust <- Mclust(precip)
plot(precipMclust)

faithfulMclust <- Mclust(faithful)
plot(faithfulMclust)

irisMclust <- Mclust(iris[,-5])
plot(irisMclust)
}
\keyword{cluster}
% docclass is function
