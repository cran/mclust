\name{emclust}
\alias{emclust}
\title{
BIC from hierarchical clustering followed by EM for several parameterized Gaussian mixture models.
}
\description{
Bayesian Information Criterion for various models and numbers of clusters
computed from hierarchical clustering followed by EM for several
parameterizations of Gaussian mixture models possibly with Poisson noise.
}
\usage{
emclust(data, nclus, modelid, k, equal=F, noise, Vinv)


}
\arguments{
\item{data}{
matrix of observations.
}
\item{nclus}{
An integer vector specifying the numbers of clusters for which the BIC is to
be calculated. Default: 1:9 without noise; 0:9 with noise.
}
\item{modelid}{
A vector of character strings indicating the models to be fitted.
The allowed values or \code{modelid} and their interpretation are as follows:
\code{"EI"} : uniform spherical, \code{"VI"} : spherical, \code{"EEE"} : uniform variance,
\code{"VVV"} : unconstrained variance, \code{"EEV"} : uniform shape and volume,
\code{"VEV"} : uniform shape.
The default is to fit all of the models.
}
\item{k}{
If \code{k} is specified, the hierarchical clustering phase will use a sample of
size \code{k} of the data in the initial hierarchical clustering phase. The
default is to use the entire data set.
}
\item{equal}{
Logical variable indicating whether or not the mixing proportions are
equal in the model. The default is to assume they are unequal.
}
\item{noise}{
A logical vector of length equal to the number of observations in the data,
whose elements indicate an initial estimate of noise (indicated by \code{T}) in
the data. By default, \code{emclust} fits Gaussian mixture models in which it is 
assumed there is no noise. If \code{noise} is specified, \code{emclust} will fit a 
Gaussian mixture with a Poisson term for noise in the EM phase.
}
\item{Vinv}{
An estimate of the inverse hypervolume of the data region (needed only if
\code{noise} is specified). Default : determined by function \code{hypvol}
}}
\value{
Bayesian Information Criterion for the six mixture models and specified
numbers of clusters. Auxiliary information returned as attributes.
}
\section{NOTE}{
The hierarchical clustering phase uses the unconstrained model.
The reciprocal condition estimate returned as an attribute ranges in value
between 0 and 1. The closer this estimate is to zero, the more likely it is
that the corresponding EM result (and BIC) are contaminated by roundoff error.
}
\references{
C. Fraley and A. E. Raftery, \emph{How many clusters? Which clustering method?
Answers via model-based cluster analysis. }Technical Report No. 329,
Dept. of Statistics, U. of Washington (February 1998).


R. Kass and A. E. Raftery, Bayes Factors. \emph{Journal of the American 
Statistical Association}90:773-795 (1995).
}
\seealso{
\code{\link{summary.emclust}}, \code{\link{emclust1}}, \code{\link{mhtree}}, \code{\link{me}}
}
\examples{
data(iris)
bicvals _ emclust(iris[,1:4], nclus=1:3, modelid=c("VVV","EEV","VEV"))


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


