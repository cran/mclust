\name{mclustBootstrapLRT}
\alias{mclustBootstrapLRT}
\alias{print.mclustBootstrapLRT}
\alias{plot.mclustBootstrapLRT}

\title{Bootstrap Likelihood Ratio Test for the Number of Mixture Components}

\description{Perform the likelihood ratio test (LRT) for assessing the number of mixture components in a specific finite mixture model parameterisation. The observed significance is approximated by using the (parametric) bootstrap for the likelihood ratio test statistic (LRTS).}

\usage{
mclustBootstrapLRT(data, modelName = NULL, nboot = 999, 
                   level = 0.05, maxG = NULL, verbose = TRUE, \dots)
                   
\method{print}{mclustBootstrapLRT}(x, \dots)

\method{plot}{mclustBootstrapLRT}(x, G = 1, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", 
    col = "forestgreen", lwd = 2, lty = 3, main = NULL, \dots) 
}

\arguments{
  \item{data}{A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables.}
  \item{modelName}{A character string indicating the mixture model to be fitted. 
  The help file for \code{\link{mclustModelNames}} describes the available models.}
  \item{nboot}{The number of bootstrap replications to use (by default 999).}
  \item{level}{The significance level to be used to terminate the sequential bootstrap procedure.}
  \item{maxG}{The maximum number of mixture components \eqn{G} to test. If not provided
  the procedure is stopped when a test is not significant at the specified \code{level}.}
  \item{verbose}{A logical, if \code{TRUE} and the session is interactive a text progress bar is displayed during the bootstrap procedure.}
  \item{\dots}{Further arguments passed to or from other methods. In particular, see the optional arguments in  \code{\link{mclustBIC}}.}
  
  \item{x}{An \code{'mclustBootstrapLRT'} object.}
  \item{G}{A value specifying the number of components for which to plot the 
           bootstrap distribution.}
  \item{hist.col}{The colour to be used to fill the bars of the histogram.}
  \item{hist.border}{The color of the border around the bars of the histogram.}
  \item{breaks}{See the argument in function \code{\link[graphics]{hist}}.}
  \item{col, lwd, lty}{The color, line width and line type to be used to represent the observed LRT statistic.}
  \item{main}{The title for the graph.}
}

\details{The implemented algorithm for computing the LRT observed significance using the bootstrap is the following.
Let \eqn{G_0} be the number of mixture components under the null hypothesis versus \eqn{G_1 = G_0+1} under the alternative. Bootstrap samples are drawn by simulating data under the null hypothesis. Then, the p-value may be approximated using eq. (13) on McLachlan and Rathnayake (2014). Equivalently, using the notation of Davison and Hinkley (1997) it may be computed as
\deqn{\textnormal{p-value} = \frac{1 + \#\{LRT^*_b \ge LRTS_{obs}\}}{B+1}}{%
      p-value = (1 + #{LRTS*_b \ge LRT_obs}) / (B+1)}
where \cr
\eqn{B} = number of bootstrap samples \cr
\eqn{LRT_{obs}}{LRT_obs} = LRTS computed on the observed data\cr
\eqn{LRT^*_b}{LRT*_b} = LRTS computed on the \eqn{b}th bootstrap sample.
}

\value{An object of class \code{'mclustBootstrapLRT'} with the following components:

  \item{G}{A vector of number of components tested under the null hypothesis.} 
  \item{modelName}{A character string specifying the mixture model as provided 
                   in the function call (see above).}
  \item{obs}{The observed values of the LRTS.}
  \item{boot}{A matrix of dimension \code{nboot} x the number of components tested 
              containing the bootstrap values of LRTS.}
  \item{p.value}{A vector of p-values.}
}

\references{
Davison, A. and Hinkley, D. (1997) \emph{Bootstrap Methods and Their Applications}. Cambridge University Press.

McLachlan G.J. (1987) On bootstrapping the likelihood ratio test statistic for the number of components in a normal mixture. \emph{Applied Statistics}, 36, 318-324.

McLachlan, G.J. and Peel, D. (2000) \emph{Finite Mixture Models}. Wiley.

McLachlan, G.J. and Rathnayake, S. (2014) On the number of components in a Gaussian mixture model. \emph{Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery}, 4(5), pp. 341-355.
}

\seealso{\code{\link{mclustBIC}}, \code{\link{mclustICL}}, \code{\link{Mclust}}}

\examples{
\dontrun{
data(faithful)
faithful.boot = mclustBootstrapLRT(faithful, model = "VVV")
faithful.boot
plot(faithful.boot, G = 1)
plot(faithful.boot, G = 2)
}
}

\keyword{htest}
\keyword{cluster}
