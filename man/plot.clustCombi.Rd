\name{plot.clustCombi}
\alias{plot.clustCombi}
\title{
Plot Combined Clusterings Results
}
\description{
Plot combined clusterings results: classifications corresponding to \code{Mclust}/BIC and to the hierarchically combined classes, and "entropy plots" to help to select a number of classes.
}
\usage{
\method{plot}{clustCombi}(x, data = NULL, what = c("classification", "entropy"), reg = 2, ...)
}
\arguments{
  \item{x}{
	Output from clustCombi.
}
  \item{data}{
	The data used to produce Output.
}
  \item{what}{
	Choose one or more of: "classification", "entropy".
}
  \item{reg}{
	The number of parts of the piecewise linear regression for the entropy plots. Choose one or more of : 2 (for 1 change-point), 3 (for 2 change-points).
}
  \item{\dots}{
	Other arguments to be passed to combiPlot, entPlot, or to the \code{Mclust} called plot functions (please see the corresponding documentations).
}
}
\value{
Classifications are plotted with combiPlot, which relies on the \code{Mclust} plot functions. 
Entropy plots (plotted with entPlot) may help to select a number of classes: please see the article cited in the references.
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
\dontrun{
data(Baudry_etal_2010_JCGS_examples)

## 1D Example 
output <- clustCombi(Test1D, G=1:15)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output, Test1D) 

## 2D Example 
output <- clustCombi(ex4.1) 

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output, ex4.1) 

## 3D Example 
output <- clustCombi(ex4.4.2)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output, ex4.4.2) 
}

}
\keyword{ cluster }
