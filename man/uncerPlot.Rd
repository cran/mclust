\name{uncerPlot}
\alias{uncerPlot}
\title{
  Uncertainty Plot for Model-Based Clustering
}
\description{
  Displays the uncertainty in converting a conditional probablility from EM
  to a classification in model-based clustering.
}
\usage{
uncerPlot(z, truth, \dots)
}
\arguments{
  \item{z}{
    A matrix whose \emph{[i,k]}th entry is the
    conditional probability of the ith observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{truth}{
    A numeric or character vector giving the true classification of the data. 
  }
  \item{\dots }{
    Provided to allow lists with elements other than the arguments can
    be passed in indirect or list calls with \code{do.call}.
  }
}
\value{
  A plot of the uncertainty profile of the data,
  with uncertainties in increasing order of magnitude.
  If \code{truth} is supplied and the number of
  classes is the same as the number of columns of 
  \code{z}, the uncertainty
  of the misclassified data is marked by vertical lines on the plot.
}
\details{
  When \code{truth} is provided and the number of classes is compatible
  with \code{z}, the function \code{compareClass} is used to to find best
  correspondence between classes in \code{truth} and \code{z}.
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
  \code{\link{mclustBIC}},
  \code{\link{em}},
  \code{\link{me}},
  \code{\link{mapClass}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
irisModel3 <-  mclustModel(iris[,-5], irisBIC, G = 3)

uncerPlot(z = irisModel3$z)
 
uncerPlot(z = irisModel3$z, truth = iris[,5])
}
\keyword{cluster}
% docclass is function
