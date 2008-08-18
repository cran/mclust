\name{plot.Mclust}
\alias{plot.Mclust}
\title{
  Plot Model-Based Clustering Results
}
\description{
  Plot model-based clustering results: BIC, classification, uncertainty and 
  (for one- and two-dimensional data) density.
}
\usage{
\method{plot}{Mclust}(x, data = NULL, what = c("BIC", "classification", 
        "uncertainty", "density"), dimens = c(1,2), xlab = NULL, ylim = NULL,
            legendArgs = list(x = "bottomright", ncol = 2, cex = 1), 
            identify = TRUE, \dots)
}
\arguments{
  \item{x}{
    Output from \code{Mclust}.
  }
  \item{data}{
    The data used to produce \code{x}.
  }
 \item{what}{
    Choose one or more of: \code{"BIC"}, \code{"classification"},   
    \code{"uncertainty"}. If the data dimension is less than 3,
    \code{"density"} can also be chosen. 
  }
\item{dimens}{
    A vector of length 2 giving the integer dimensions of the
    desired coordinate projections for multivariate data. The default is
    \code{c(1,2)}, in which the first dimension is plotted against the second.
  }
 \item{xlab}{
    Optional label for the horizontal axis of the BIC plot.
  }
 \item{ylim}{
    Optiomnal limits for the vertical axis of the BIC plot.
  }
 \item{legendArgs}{
    Arguments to pass to the \code{legend} function. Set to \code{NULL}
    for no legend.
  }
 \item{identify}{
    A logical variable indicating whether or not to add a title to the plot
    identifying the dimensions used.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\value{
  Model-based clustering plots: BIC values used for choosing the number of 
  clusters. For data in more than two dimensions, a pairs plot of the showing 
  the classification, a coordinate projections of the data showing location of 
  the mixture components, classification, and uncertainty. For one- and two- 
  dimensional  data, plots showing location of the mixture components,
  classification, uncertainty, and density.  
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

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}

\seealso{
  \code{\link{Mclust}},
  \code{\link{mclust1Dplot}},
  \code{\link{mclust2Dplot}},
  \code{\link{surfacePlot}},
  \code{\link{coordProj}},
  \code{\link{randProj}}
}
\examples{
\dontrun{
plot(Mclust(precip),precip)

plot(Mclust(faithful),faithful)

plot(Mclust(iris[,-5]),iris[,-5])
}
}
\keyword{cluster}
% docclass is function
