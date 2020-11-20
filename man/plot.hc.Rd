\name{plot.hc}
\alias{plot.hc}

\title{Dendrograms for Model-based Agglomerative Hierarchical Clustering}

\description{
   Display two types for dendrograms for model-based hierarchical clustering
   objects.
}

\usage{
\method{plot}{hc}(x, what=c("loglik","merge"), maxG=NULL, labels=FALSE, hang=0, \dots)
}

\arguments{
  \item{x}{
     An object of class \code{'hc'}.
  }
  \item{what}{
    A character string indicating the type of dendrogram to be displayed.\cr
    Possible options are:
    \describe{
    \item{\code{"loglik"}}{Distances between dendrogram levels are based on
          the classification likelihood.}
    \item{\code{"merge"}}{Distances between dendrogram levels are uniform,
                          so that levels correspond to the number of clusters.}
    }
  }
  \item{maxG}{
    The maximum number of clusters for the dendrogram.
    For \code{what = "merge"}, the default is the
    number of clusters in the initial partition.
    For \code{what = "loglik"}, the default is the minimnum of the
    maximum number of clusters for which the classification loglikelihood
    an be computed in most cases, and the maximum number of clusters for
    which the classification likelihood increases with increasing numbers of
    clusters.
  }
  \item{labels}{
    A logical variable indicating whether or not to display leaf (observation)
    labels for the dendrogram (row names of the data). These are likely to be 
    useful only if the number of observations in fairly small, since otherwise
    the labels will be too crowded to read. 
    The default is not to display the leaf labels.
  }
  \item{hang}{
    For \code{hclust} objects, this argument is the fraction of the plot 
    height by which labels should hang below the rest of the plot. A negative 
    value will cause the labels to hang down from 0.
    Because model-based hierarchical clustering does not share all of the
    properties of \code{hclust}, the \code{hang} argment won't work in
    many instances.
  }
  \item{\dots}{
    Additional plotting arguments.
  }
}
\value{
  A dendrogram is drawn, with distances based on either the classification
  likelihood or the merge level (number of clusters). 
}
\details{
   The plotting input does not share all of the properties of \code{hclust}
   objects, hence not all plotting arguments associated with \code{hclust}
   can be expected to work here.
}
\note{
  If \code{modelName = "E"} (univariate with equal variances) or
  \code{modelName = "EII"} (multivariate with equal spherical
  covariances), then the underlying model is the same as for
  Ward's method for hierarchical clustering.
}
\references{
  J. D. Banfield and A. E. Raftery (1993).
  Model-based Gaussian and non-Gaussian Clustering.
  \emph{Biometrics 49:803-821}. 
  
  C. Fraley (1998).
  Algorithms for model-based Gaussian hierarchical clustering.
  \emph{SIAM Journal on Scientific Computing 20:270-281}. 
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
}
\seealso{
  \code{\link{hc}}
}
\examples{
data(EuroUnemployment)
hcTree <- hc(modelName = "VVV", data = EuroUnemployment)
plot(hcTree, what = "loglik")
plot(hcTree, what = "loglik", labels = TRUE)
plot(hcTree, what = "loglik", maxG = 5, labels = TRUE)
plot(hcTree, what = "merge")
plot(hcTree, what = "merge", labels = TRUE)
plot(hcTree, what = "merge", labels = TRUE, hang = 0.1)
plot(hcTree, what = "merge", labels = TRUE, hang = -1)
plot(hcTree, what = "merge", labels = TRUE, maxG = 5)
}
\keyword{cluster}
