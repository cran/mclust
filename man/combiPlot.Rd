\name{combiPlot}
\alias{combiPlot}
\title{
Plot Classifications Corresponding to Successive Combined Solutions
}
\description{
Plot classifications corresponding to successive combined solutions. 
}
\usage{
combiPlot(data, z, combiM, ...)
}
\arguments{
  \item{data}{
	The data.
}
  \item{z}{
	A matrix whose [i,k]th entry is the probability that observation i in the data belongs to the kth class, for the initial solution (ie before any combining). Typically, the one returned by \code{Mclust}/BIC.
}
  \item{combiM}{
	A "combining matrix" (as provided by \code{\link{clustCombi}}), ie a matrix whose kth row contains only zeros, but in columns corresponding to the labels of the classes in the initial solution to be merged together to get the combined solution.
}
  \item{\dots}{
	Other arguments to be passed to the \code{\link{Mclust}} plot functions.
}
}
\value{
Plot the classifications obtained by MAP from the matrix \code{t(combiM \%*\% t(z))}, which is the matrix whose [i,k]th entry is the probability that observation i in the data belongs to the kth class, according to the combined solution obtained by merging (according to \code{combiM}) the initial solution described by \code{z}.
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{clustCombi}}, \code{\link{combMat}}, \code{\link{clustCombi}}
}
\examples{

data(Baudry_etal_2010_JCGS_examples)
MclustOutput <- Mclust(ex4.1) 

MclustOutput$G # Mclust/BIC selected 6 classes

par(mfrow=c(2,2))

combiM0 = diag(6) # is the identity matrix
# no merging: plot the initial solution, given by z
combiPlot(ex4.1, MclustOutput$z, combiM0, cex = 3) 
title("No combining")

combiM1 = combMat(6, 1, 2) # let's merge classes labeled 1 and 2
combiM1
combiPlot(ex4.1, MclustOutput$z, combiM1)
title("Combine 1 and 2")

# let's merge classes labeled 1 and 2, and then components labeled (in this 
# new 5-classes combined solution...) 1 and 2
combiM2 = combMat(5, 1, 2) \%*\% combMat(6, 1, 2) 
combiM2 
combiPlot(ex4.1, MclustOutput$z, combiM2)
title("Combine 1, 2 and then 1 and 2 again")

plot(rep(1,6), 1:6, col = mclust.options()$classPlotColors, 
     pch = mclust.options()$classPlotSymbols, 
     xlab = "", ylab = "Class label", xaxt = "n")
title("legend")
}
\keyword{ cluster }
