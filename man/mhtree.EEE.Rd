\name{mhtree.EEE}
\alias{mhtree.EEE}
\title{
Classification tree for hierarchical clustering for Gaussian models with constant variance.
}
\section{DESCRIPTION}{
Computes a classification tree for agglomerative hierarchical clustering using
a Gaussian model in which all clusters have the same volume, shape, and 
orientation.
}
\usage{
mhtree.EEE(data, partition, min.clusters = 1)


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
an object of class \code{"mhclust"}, which consists of a classification tree with
the following attributes:

\item{call}{
a copy of the call to \code{mhclust.EEE}.
}
\item{determinant}{
value related to the determinant of the sum of sample crossproduct matrices,
which is the relevant criterion for the model, computed during the course of 
hierarchical clustering.
}
\item{dimensions}{
the data dimensions.
}
\item{initial.partition}{
the partition at which agglomerative hierarchical clustering is initiated.
}
\item{trace}{
the trace of the sum of sample crossproduct matrices produced during the
course of hierarchical clustering.
}}
\section{NOTES}{
The constant variance option is one of the slowest options for \code{mhclust},
because the model does not admit a fast hierarchical clustering algorithm.
At the same time, it is one of the more space-sefficient options.
}
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.


H. P. Friedman and J. Rubin, On some invariant criteria for grouping data,
\emph{Journal of the American Statistical Association, }62:1159-1178 (1967).
}
\seealso{
\code{\link{mhtree}}, \code{\link{mhclass}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)
mhtree.EEE(iris[,1:4])


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


