\name{clustCombiOptim}
\alias{clustCombiOptim}

\title{Optimal number of clusters obtained by combining mixture components}

\description{
Return the optimal number of clusters by combining mixture components based on the entropy method discussed in the reference given below.
}

\usage{
clustCombiOptim(object, reg = 2, plot = FALSE, \dots)
}
\arguments{
  \item{object}{
  An object of class \code{'clustCombi'} resulting from a call to \code{\link{clustCombi}}.
}
  \item{reg}{
  The number of parts of the piecewise linear regression for the entropy plots. 
  Choose 2 for a two-segment piecewise linear regression model (i.e. 1 change-point), and 3 for a three-segment piecewise linear regression model (i.e. 3 change-points).
}
  \item{plot}{
   Logical, if \code{TRUE} an entropy plot is also produced.
}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function returns a list with the following components:
\item{numClusters.combi}{The estimated number of clusters.}
\item{z.combi}{A matrix whose \emph{[i,k]}th entry is the probability that observation \emph{i} in the data belongs to the \emph{k}th cluster.}
\item{cluster.combi}{The clustering labels.}
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{combiPlot}}, \code{\link{entPlot}}, \code{\link{clustCombi}}
}
\examples{
data(Baudry_etal_2010_JCGS_examples)
output <- clustCombi(data = ex4.1) 
combiOptim <- clustCombiOptim(output)
str(combiOptim)

# plot optimal clustering with alpha color transparency proportional to uncertainty
zmax <- apply(combiOptim$z.combi, 1, max)
col <- mclust.options("classPlotColors")[combiOptim$cluster.combi]
vadjustcolor <- Vectorize(adjustcolor)
alphacol = (zmax - 1/combiOptim$numClusters.combi)/(1-1/combiOptim$numClusters.combi)
col <- vadjustcolor(col, alpha.f = alphacol)
plot(ex4.1, col = col, pch = mclust.options("classPlotSymbols")[combiOptim$cluster.combi])
}

\keyword{ cluster }
