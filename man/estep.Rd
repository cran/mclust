\name{estep}
\alias{estep}
\title{
E-step for parameterized MVN mixture models
}
\description{
E-step for estimating conditional probabilities from parameter estimates in an
MVN mixture model having possibly one Poisson noise term.
}
\usage{
estep(data, modelid, mu, \dots)


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
\item{mu}{
matrix whose columns are the Gaussian group means.
}
\item{...}{
additional arguments, as follows:
}
\item{sigmasq}{
}
\item{sigma}{
group variances (\code{sigmasq} - spherical models) or covariances (\code{sigma} -
elliposidal models)
}
\item{prob}{
mixing proportions (probabilities) for each group. If \code{prob} is missing, 
the number of groups is assumed to be the number of columns in \code{mu} (no
noise). A Poisson noise term will appear in the conditional probabilities if
\code{length(prob)} is equal to \code{ncol(mu)+1}.
}
\item{eps}{
Tolerance for determining singularity in the covariance matrix. The precise
definition of \code{eps} varies the parameterization, each of which has a default.
}
\item{Vinv}{
An estimate of the inverse hypervolume of the data region (needed only if
\code{prob} indicates a noise term). Default : determined by function \code{hypvol}
}}
\value{
the conditional probablilities corresponding to the parameter estimates.
The loglikelihood is returned as an attribute.
}
\section{NOTE}{
The reciprocal condition estimate returned as an attribute ranges in value
between 0 and 1. The closer this estimate is to zero, the more likely it is
that the corresponding EM result (and BIC) are contaminated by roundoff error.
}
\references{
G. Celeux and G. Govaert, Gaussian parsimonious clustering models,
\emph{Pattern Recognition, }28:781-793 (1995).


A. P. Dempster, N. M. Laird and D. B. Rubin, Maximum Likelihood from
Incomplete Data via the EM Algorithm, \emph{Journal of the Royal Statistical
Society, Series B, }39:1-22 (1977).


G. J. MacLachlan and K. E. Basford, The EM Algorithm and Extensions, Wiley
(1997).
}
\seealso{
\code{\link{me}}, \code{\link{mstep}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid="VI"), 3)
z <- me( iris[,1:4], ctoz(cl), modelid = "VI")
pars <- mstep( iris[,1:4], modelid = "VI", z)
estep(iris[,1:4], modelid = "VI", pars$mu, pars$sigma, pars$prob)


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


