\name{mstep.EI}
\alias{mstep.EI}
\title{
M-step for spherical, constant-volume MVN mixture models
}
\description{
M-step for estimating parameters given conditional probabilities in an MVN 
mixture model having equal, spherical variances and possibly one Poisson noise
term.
}
\usage{
mstep.EI(data, z, eps, equal = F, noise = F, Vinv)


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
Lower bound on the estimated value of sigma-squared.
Default : \code{.Machine$double.eps}
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
\code{mstep}, \code{me.EI}, \code{estep.EI}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid = "EI"),3)
z <- me.EI( iris[,1:4], ctoz(cl))
Mstep <- mstep.EI(iris[,1:4], z)
estep.EI( iris[,1:4], Mstep$mu, Mstep$sigma, Mstep$prob)


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


