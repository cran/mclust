\name{plot.clustCombi}
\alias{plot.clustCombi}
\title{
Plot Combined Clusterings Results
}
\description{
Plot combined clusterings results: classifications corresponding to \code{Mclust}/BIC and to the hierarchically combined classes, "entropy plots" to help to select a number of classes, and the tree structure obtained from combining mixture components.
}
\usage{
\method{plot}{clustCombi}(x, what = c("classification", "entropy", "tree"), \dots)
}
\arguments{
  \item{x}{
  Object returned by \code{\link{clustCombi}} function.
}
  \item{what}{
  Type of plot.
}
  \item{\dots}{
  Other arguments to be passed to other functions: \code{\link{combiPlot}}, \code{\link{entPlot}}, \code{\link{combiTree}}. Please see the corresponding documentations.
}
}
\value{
Classifications are plotted with \code{\link{combiPlot}}, which relies on the \code{Mclust} plot functions. 
Entropy plots are plotted with \code{\link{entPlot}} and may help to select a number of classes: please see the article cited in the references.
Tree plots are produced by \code{\link{combiTree}} and graph the tree structure implied by the clusters combining  process. 
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{combiPlot}}, \code{\link{entPlot}}, \code{\link{combiTree}}, \code{\link{clustCombi}}.
}
\examples{
\dontrun{
data(Baudry_etal_2010_JCGS_examples)

## 1D Example 
output <- clustCombi(data = Test1D, G=1:15)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output) 

## 2D Example 
output <- clustCombi(data = ex4.1) 

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output) 

## 3D Example 
output <- clustCombi(data = ex4.4.2)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output)
}

}
\keyword{ cluster }
