\name{mstep.VEV}
\alias{mstep.VEV}
\title{
M-step for constant shape, constant volume MVN mixture models
}
\description{
M-step for estimating parameters given conditional probabilities in an MVN 
mixture model having constant shape, constant volume and possibly one Poisson 
noise term.
}
\usage{
mstep.VEV(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)


}
\arguments{
\item{data}{
matrix of observations.
}
\item{z}{
matrix of conditional probabilities. \code{z} should have a row for each observation
in \code{data}, and a column for each component of the mixture.
}
\item{eps}{
A 2-vector specifying lower bounds on the pth root of the volume of the 
ellipsoids defining the clusters, where p is the data dimension, and on the 
reciprocal condition number for the estimated shape of the covariance 
estimates. Default: \code{c(.Machine$double.eps, .Machine$double.eps)}
}
\item{tol}{
The iteration for volume/shape estimates is terminated if their relative
error is less than \code{tol}.
}
\item{itmax}{
The iteration for volume/shape estimates is terminated if the number of 
iterations exceeds \code{itmax}. Default: \code{Inf} (termination is determined by
\code{tol}).
}
\item{equal}{
Logical variable indicating whether or not to assume equal proportions in the
mixture. Default : \code{F}.
}
\item{noise}{
Logical variable indicating whether or not to include a Poisson noise term in
the model. Default : \code{F}.
}
\item{Vinv}{
An estimate of the inverse hypervolume of the data region (needed only if
\code{noise = T}). Default : determined by function \code{hypvol}
}}
\value{
A list whose components are the parameter estimates corresponding to \code{z}:

\item{mu}{
matrix whose columns are the Gaussian group means.
}
\item{sigma}{
group variance matrix.
}
\item{prob}{
probabilities (mixing proportions) for each group (present only when
\code{equal = T}).
The loglikelihood and reciprocal condition estimate are returned as attributes.
}}
\references{
G. Celeux and G. Govaert, Gaussian parsimonious clustering models,
\emph{Pattern Recognition, }28:781-793 (1995).


A. P. Dempster, N. M. Laird and D. B. Rubin, Maximum Likelihood from
Incomplete Data via the EM Algorithm, \emph{Journal of the Royal Statistical
Society, Series B, }39:1-22 (1977).


G. J. MacLachlan and K. E. Basford, The EM Algorithm and Extensions, Wiley,
(1997).
}
\seealso{
\code{\link{mstep}}, \code{\link{me.VEV}}, \code{\link{estep.XEV}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4]),3)
z <- me.VEV( iris[,1:4], ctoz(cl))
mstep.VEV(iris[,1:4], z)


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


