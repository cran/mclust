\name{mhtree.EI}
\alias{mhtree.EI}
\title{
Classification tree for hierarchical clustering for Gaussian models with uniform diagonal variance.
}
\description{
Computes a classification tree for agglomerative hierarchical clustering using
a Gaussian model in which clusters are spherical and of equal volume (Wards' 
method).
}
\usage{
mhtree.EI(data, partition, min.clusters = 1)


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
}}
\value{
an object of class \code{"mhtree"}, which consists of a classification tree with
the following attributes:

\item{call}{
a copy of the call to \code{mhtree.EI}.
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
Clustering, \emph{Biometrics, }49:803-821 (September 1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.


J. H. Ward, Hierarchical groupings to optimize and objective function,
\emph{Journal of the American Statistical Association, }58:234-244 (1963).
}
\seealso{
\code{\link{mhtree}}, \code{\link{mhclass}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)
mhtree.EI(iris[,1:4])


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


