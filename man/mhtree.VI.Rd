\name{mhtree.VI}
\alias{mhtree.VI}
\title{
Classification tree for hierarchical clustering for Gaussian models with diagonal variance.
}
\description{
Computes a classification tree for agglomerative hierarchical clustering using 
a Gaussian model in
which clusters are spherical and may vary in volume volume.
}
\usage{
mhtree.VI(data, partition, min.clusters = 1, alpha = 1)


}
\arguments{
\item{data}{
matrix of observations.
}
\item{partition}{
initial classification of the data. The default puts every observation in a
singleton cluster.
}
\item{min.clusters}{
minimum number of clusters desired. The default is to carry out agglomerative
hierarchical clustering until termination, that is, until all observations
belong to a single group. The default value is 1.
}
\item{alpha}{
The value \code{alpha} times the trace of the sample crossproduct matrix of all
the observations divided by the product of the data dimensions, is used for
the purpose of initalization.
The default value is 1.
}}
\value{
an object of class \code{"mhtree"}, which consists of a classification tree with
the following attributes:

\item{call}{
a copy of the call to \code{mhtree.VI}.
}
\item{change}{
value of the optimal change in likelihood at each stage.
}
\item{dimensions}{
the data dimensions.
}
\item{initial.partition}{
the partition at which agglomerative hierarchical clustering is initiated.
}}
\section{NOTES}{
The value \code{alpha} is needed because the criterion is not defined for singleton
clusters or clusters consisting only of points that coincide.
}
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.
}
\seealso{
\code{\link{mhtree}}, \code{\link{mhclass}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)
mhtree.VI(iris[,1:4])


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


