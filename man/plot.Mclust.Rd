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
plot.Mclust(x, data, dimens = c(1, 2), scale = FALSE, \dots)
}
\arguments{
  \item{x}{
    Output from \code{Mclust}.
  }
  \item{data}{
    The data used to produce \code{x}.
  }
  \item{dimens}{
    An integer vector of length two specifying the dimensions for
    coordinate projections if the data is more than two-dimensional. The
    default is \code{c(1,2)} (the first two dimesions). 
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and thus preserve
    the shape of the distribution. Default: \code{scale=FALSE}  
  }
  \item{\dots}{
    Further arguments to the lower level plotting functions.
  }
  
}
\value{
  Plots selected via a menu including the following options: BIC values
  used for choosing the number of clusters For data in more than two
  dimensions, a pairs plot of the showing the classification, coordinate
  projections of the data, showing location of the mixture components,
  classification, and/or uncertainty. For one- and two- dimensional
  data, plots showing location of the mixture components,
  classification, uncertainty, and or density.  
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}.
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington.
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{Mclust}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisMclust <- Mclust(irisMatrix)

\dontrun{plot(irisMclust,irisMatrix)}
}
\keyword{cluster}
% docclass is function
