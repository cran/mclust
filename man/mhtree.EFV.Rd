\name{mhtree.EFV}
\alias{mhtree.EFV}
\title{
Classification tree for hierarchical clustering for Gaussian models with equal volume and fixed shape.
}
\description{
Computes a classicification tree for agglomerative hierarchical clustering 
using a Gaussian model in which clusters have a fixed shape and are of equal 
volume.
}
\usage{
mhtree.EFV(data, partition, min.clusters = 1, shape)


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
\item{shape}{
positive vector of length equal to the dimension of the observations.
The values of \code{shape} should be proportional to the squares of the lengths
of the semi-major axes of the corresponding ellipsoid, and may be given in
any order. Note that \code{shape} is a required argument.
}}
\value{
an object of class \code{"mhtree"}, which consists of a classification tree with
the following attributes:

\item{call}{
a copy of the call to \code{mhtree.EFV}.
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
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.


F. Murtagh and A. E. Raftery, Fitting straight lines to point patterns,
\emph{Pattern Recognition, }17:479-483 (1984).
}
\seealso{
\code{\link{mhtree}}, \code{\link{mhclass}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)
shape <- c(1,1/2,1/3)
mhtree.EFV(iris[,1:3], shape = shape)


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


