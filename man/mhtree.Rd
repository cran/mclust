\name{mhtree}
\alias{print.mhtree}
\title{
Classification Tree for Model-based Gaussian hierarchical clustering.
}
\description{
Determines a classification tree for agglomerative hierarchical
clustering using criteria based on parameterizations of Gaussian mixture
models that reflect the underlying geometry of the resulting clusters.
}
\usage{
mhtree(data, modelid, partition, min.clusters = 1, verbose = F, \dots)
print.mhtree(tree.object, \dots)

}
\arguments{
\item{data}{
matrix of observations.
}
\item{modelid}{
An integer specifying a parameterization of the MVN covariance matrix defined 
by volume, shape and orientation charactertistics of the underlying clusters. 
The allowed values or \code{model} and their interpretation are as follows:
\code{"EI"} : uniform spherical, \code{"VI"} : spherical, \code{"EEE"} : uniform variance, 
\code{"VVV"} : unconstrained  variance, \code{"EFV"} : fixed (user supplied)
uniform volume, \code{"VFV"} : fixed (user supplied) shape.
}
\item{partition}{
initial classification of the data. The default puts every observation in a
singleton cluster.
}
\item{min.clusters}{
minimum number of clusters desired. The default is to carry out agglomerative
hierarchical clustering until termination, that is, until all observations
belong to a single group.
}
\item{verbose}{
A logical variable specifying printing of the model type when set to \code{T}.
}
\item{...}{
Allows users to specify the required \code{shape} parameter for the two fixed
shape models \code{"EFV"} and \code{"VFV"}, and to change default parameters that are 
used in the algorithms underlying some models.
}}
\value{
an object of class \code{"mhclust"}, which consists of a classification tree with
attributes giving other information relating to the clustering process.
}
\section{NOTES}{
Only the six models illustrated in the example below are supported at present.
These correspond to the models discussed in the Banfield and Raftery 
reference.
It may be desirable to transform the data in some way before attempting to 
partition it into clusters.
Different permuations of the data may produce different classifications,
because \code{mhclust} resolves ties in a way that is dependent on the order of the
observations, and because values of criterion that are close may change enough
to affect the choice of merge pairs in a given stage.
}
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (September 1993).


C. Fraley, \emph{Algorithms for Model-based Gaussian Hierarchical Clustering,
}Technical Report No. 311, Department of Statistics, University of 
Washington (October 1996), to appear in
\emph{SIAM Journal on Scientific Computing}.
}
\seealso{
\code{\link{mhclass}}, \code{\link{loglik}}, \code{\link{awe}}, \code{\link{partuniq}}
}
\examples{
data(iris)

# Ellipsoidal, equal volume, shape and orientation
mhtree(iris[,1:4], modelid = "EEE")

# Spherical, equal volume, fixed shape, variable orientation
shape <- c(1,1/2,1/3,1/4)
mhtree(iris[,1:4], modelid = "EFV", shape=shape)

# Spherical, equal volume (Ward's method).
mhtree(iris[,1:4], modelid = "EI")

# Ellipsoidal, equal volume, constant shape, variable orientation
mhtree(iris[,1:4], modelid = "VFV", shape=shape)

# Spherical, variable volume
mhtree(iris[,1:4], modelid = "VI")

# Unconstrained (default).
mhtree(iris[,1:4], modelid = "VVV")

}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


