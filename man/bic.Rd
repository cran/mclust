\name{bic}
\alias{bic.EI}
\alias{bic.VI}
\alias{bic.EEE}
\alias{bic.EEV}
\alias{bic.VEV}
\alias{bic.VVV}
\title{
BIC for parameterized MVN mixture models
}
\description{
Bayesian Information Criterion for MVN mixture models with possibly one 
Poisson noise term.
}
\usage{
bic(data, modelid, z, eps, tol, itmax, equal = F, noise = F, Vinv)


}
\arguments{
\item{data}{
matrix of observations.
}
\item{modelid}{
An integer specifying a parameterization of the MVN covariance matrix defined 
by volume, shape and orientation charactertistics of the underlying clusters. 
The allowed values for \code{modelid} and their interpretation are as follows:
\code{"EI"} : uniform spherical, \code{"VI"} : spherical, \code{"EEE"} : uniform variance, 
\code{"VVV"} : unconstrained variance, \code{"EEV"} : uniform shape and volume, 
\code{"VEV"} : uniform shape.
}
\item{...}{
other arguments, including a quantity \code{eps} for determining singularity
in the covariance, and the following:
}
\item{z}{
matrix of conditional probabilities. \code{z} should have a row for each observation
in \code{data}, and a column for each component of the mixture. If \code{z} is missing,
a single cluster is assumed (all noise if \code{noise = T}).
}
\item{eps}{
Tolerance for determining singularity in the covariance matrix. The precise 
definition of \code{eps} varies the parameterization, each of which has a default.
}
\item{equal}{
Logical variable indicating whether or not the mixing proportions are
equal in the model. The default is to assume they are unequal.
}
\item{noise}{
Logical variable indicating whether or not to include a Poisson noise term in
the model. Default : \code{F}.
}
\item{Vinv}{
An estimate of the inverse hypervolume of the data region (needed only if
\code{noise = T}). Default : determined by the function \code{hypvol}
}}
\value{
An object of class \code{"bic"} which is the Bayesian Information Criterion for the
given mixture model and given conditional probabilites. The model parameters 
and reciprocal condition estimate are returned as attributes.
}
\section{NOTE}{
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
\code{\link{me}}, \code{\link{mstep}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid = "VVV"), 3)
z <- me( iris[,1:4], ctoz(cl), modelid = "VVV")
bic(iris[,1:4], modelid = "VVV", z = z)


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


