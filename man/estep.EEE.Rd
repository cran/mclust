\name{estep.EEE}
\alias{estep.EEE}
\title{
E-step for constant-variance MVN mixture models
}
\description{
E-step for estimating conditional probabilities from parameter estimates in an
MVN mixture model having constant variance and possibly one Poisson noise term.
}
\usage{
estep.EEE(data, mu, sigma, prob, eps, Vinv)


}
\arguments{
\item{data}{
matrix of observations.
}
\item{mu}{
matrix whose columns are the Gaussian group means.
}
\item{sigma}{
group variance matrix.
}
\item{prob}{
mixing proportions (probabilities) for each group. If \code{prob} is missing, 
the number of groups is assumed to be the number of columns in \code{mu} (no
noise). A Poisson noise term will appear in the conditional probabilities if
\code{length(prob)} is equal to \code{ncol(mu)+1}.
}
\item{eps}{
Lower bound on the reciprocal condition estimate for the covariance. 
Default : \code{.Machine$double.eps}.
}
\item{Vinv}{
An estimate of the inverse hypervolume of the data region (needed only if
\code{prob} indicates a noise term). Default : determined by function \code{hypvol}
}}
\value{
the conditional probablilities corresponding to the parameter estimates.
The loglikelihood is returned as an attribute.
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
\code{\link{estep}}, \code{\link{me.EEE}}, \code{\link{mstep.EEE}}
}
\examples{
data(iris)
cl <- mhclass(mhtree(iris[,1:4], modelid = "EEE"),3)
z <- me.EEE( iris[,1:4], ctoz(cl))
Mstep <- mstep.EEE(iris[,1:4], z)
estep.EEE( iris[,1:4], Mstep$mu, Mstep$sigma, Mstep$prob)


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


