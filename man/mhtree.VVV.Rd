\name{mhtree.VVV}
\alias{mhtree.VVV}
\title{
Classification tree for hierarchical clustering for Gaussian models with unconstrained variance.
}
\section{DESCRIPTION}{
Computes a classification tree for agglomerative hierarchical clustering using
a Gaussian model in which volume, shape, and orientation may vary among 
clusters.
}
\usage{
mhtree.VVV(data, partition, min.clusters = 1, alpha = 1, beta = 1)


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
}
\item{beta}{
The value \code{beta} times the trace of the sample crossproduct matrix for the
group divided by the number of observations in the group, is used for the
purpose of initialization.
The default value is 1.
}}
\value{
an object of class \code{"mhtree"}, which consists of a classification tree with
the following attributes:

\item{call}{
a copy of the call to \code{mhtree.VVV}.
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
The values \code{alpha} and \code{beta} are needed since the criteria are not
defined for groups for which the rank of the sample cross product matrix is
less than the dimension of the observations.
}
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.


A. J. Scott and M. J. Simons, Clustering methods based on likelihood ratio
criteria, \emph{Biometrics, }27:387-397 (1971).
}
\seealso{
\code{\link{mhtree}}, \code{\link{mhclass}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)
mhtree.VVV(iris[,1:4])


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


