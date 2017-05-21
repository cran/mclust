\name{gmmhd}
\alias{gmmhd}
\alias{print.gmmhd}
\alias{summary.gmmhd}
\alias{print.summary.gmmhd}
\alias{plot.gmmhd}
\alias{gmmhdClusterCores}
\alias{gmmhdClassify}
  
\title{Identifying Connected Components in Gaussian Finite Mixture Models for Clustering}

\description{
Starting with the density estimate obtained from a fitted Gaussian finite mixture model, cluster cores are identified from the connected components at a given density level. Once cluster cores are identified, the remaining observations are allocated to those cluster cores for which the probability of cluster membership is the highest. 
}

\usage{
gmmhd(object, 
      ngrid = min(round((log(nrow(data)))*10), nrow(data)), 
      dr = list(d = 3, lambda = 1, cumEvalues = NULL, mindir = 2),
      classify = list(G = 1:5, 
                      modelNames = mclust.options("emModelNames")[-c(8, 10)]),
      \dots)

\method{plot}{gmmhd}(x, what = c("mode", "cores", "clusters"), \dots)
}

\arguments{
  \item{object}{An object returned by \code{\link{Mclust}}.}
  \item{ngrid}{An integer specifying the number of grid points used to compute the density levels.}
  \item{dr}{A list of parameters used in the dimension reduction step.}
  \item{classify}{A list of parameters used in the classification step.}
  \item{x}{An object of class \code{'gmmhd'} as returned by the function \code{gmmhd}.}
  \item{what}{A string specifying the type of plot to be produced. See Examples section.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
Model-based clustering associates each component of a finite mixture distribution to a group or cluster.
An underlying implicit assumption is that a one-to-one correspondence exists between mixture components and clusters. However, a single Gaussian density may not be sufficient, and two or more mixture components could be needed to reasonably approximate the distribution within a homogeneous group of observations. 

This function implements the methodology proposed by Scrucca (2016) based on the identification of high density regions of the underlying density function. Starting with an estimated Gaussian finite mixture model, the corresponding density estimate is used to identify the cluster cores, i.e. those data points which form the core of the clusters. 
These cluster cores are obtained from the connected components at a given density level \eqn{c}. A mode function gives the number of connected components as the level \eqn{c} is varied. 
Once cluster cores are identified, the remaining observations are allocated to those cluster cores for which the probability of cluster membership is the highest.

The method usually improves the identification of non-Gaussian clusters compared to a fully parametric approach. Furthermore, it enables the identification of clusters which cannot be obtained by merging mixture components, and it can be straightforwardly extended to cases of higher dimensionality.
}
\value{
A list of class \code{gmmhd} with the following components:

\item{Mclust}{The input object of class \code{"Mclust"} representing an estimated Gaussian finite mixture model.}
\item{MclustDA}{An object of class \code{"MclustDA"} containing the model used for the classification step.}
\item{MclustDR}{An object of class \code{"MclustDR"} containing the dimension reduction step if performed, otherwise \code{NULL}.}
\item{x}{The data used in the algorithm. This can be the input data or a projection if a preliminary dimension reduction step is performed.}
\item{density}{The density estimated from the input Gaussian finite mixture model evaluated at the input data.}
\item{con}{A list of connected components at each step.}
\item{nc}{A vector giving the number of connected components (i.e. modes) at each step.}
\item{pn}{Vector of values over a uniform grid of proportions of length \code{ngrid}.}
\item{qn}{Vector of density quantiles corresponding to proportions \code{pn}.}
\item{pc}{Vector of empirical proportions corresponding to quantiles \code{qn}.}
\item{clusterCores}{Vector of cluster cores numerical labels; \code{NA}s indicate that an observation does not belong to any cluster core.}
\item{clusterCores}{Vector of numerical labels giving the final clustering.}
\item{numClusters}{An integer giving the number of clusters.}
}

\references{
Scrucca, L. (2016) Identifying connected components in Gaussian finite mixture models for clustering. \emph{Computational Statistics & Data Analysis}, 93, 5-17.
}

\author{
Luca Scrucca \email{luca.scrucca@unipg.it}
}

%\note{}

\seealso{\code{\link{Mclust}}}

\examples{
\dontrun{
data(faithful)
mod <- Mclust(faithful)
summary(mod)
plot(as.densityMclust(mod), faithful, what = "density", 
     points.pch = mclust.options("classPlotSymbols")[mod$classification], 
     points.col = mclust.options("classPlotColors")[mod$classification])

GMMHD <- gmmhd(mod)
summary(GMMHD)

plot(GMMHD, what = "mode")
plot(GMMHD, what = "cores")
plot(GMMHD, what = "clusters")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. Use one of  RShowDoc("KEYWORDS")
\keyword{cluster}
