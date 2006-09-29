\name{plot.mclustDAtrain}
\alias{plot.mclustDAtrain}
\title{
  Plot mclustDA training models.
}
\description{
  Plots representation of the models produced by
  \code{mclustDAtrain}. For multidimensional data,
   the plot is a coordinate projection and the ellipses shown correspond to 
   the covariance matrices.
}
\usage{
plot.mclustDAtrain(x, data, dimens=c(1,2), symbols=NULL, 
          scale = FALSE, xlim=NULL, ylim=NULL, CEX = 1, \dots)
}
\arguments{
  \item{x}{
    An object produced by a call to \code{mclustDAtrain}.
  }
  \item{data}{
    A numeric matrix or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{dimens}{
    A vector of length 2 giving the integer dimensions of the
    desired coordinate projections. The default is
    \code{c(1,2)}, in which the first
    dimension is plotted against the second.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to
    each unique class in \code{classification}. Elements in
    \code{symbols} correspond to classes in \code{classification} in
    sorted order. Default: If \emph{G} is the number of groups in the
    classification, the first \emph{G} symbols in 
     \code{.Mclust\$classPlotSymbols}, otherwise if \emph{G} is less than 
     27 then the first \emph{G} capital letters in the Roman alphabet. 
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=FALSE} 
  }
  \item{xlim, ylim}{
    Arguments specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\section{Side Effects}{
  A plot showing a two-dimensional coordinate projection of the data, together
  with the location of the  mixture components, classification, uncertainty,
  and/or classification errors.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{coordProj}},
  \code{\link{mclust1Dplot}},
  \code{\link{mclust2Dplot}},
  \code{\link{mclustOptions}}
}
\examples{
odd <- seq(from = 1, to = nrow(iris), by = 2)

irisTrain <- mclustDAtrain(data = iris[odd,-5], labels = iris[odd,5])
\dontrun{
plot(irisTrain, iris[odd,-5])
}
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
