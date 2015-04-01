\name{print.clustCombi}
\alias{print.clustCombi}
\title{
Displays Combined Clusterings Results
}
\description{
Displays the Gaussian mixture model and number of components selected by Mclust/BIC and the components proportions and means of the fitted corresponding distribution. Then, informations about the combining steps.
}
\usage{
\method{print}{clustCombi}(x, \dots) 
}
\arguments{
  \item{x}{An object of class \code{'clustCombi'}. Typically an output from the \code{clustCombi} function.}
\item{\dots}{additional arguments.}
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{clustCombi}}	
}
\examples{
\dontrun{
data(Baudry_etal_2010_JCGS_examples)

## 2D Example 
output <- clustCombi(ex4.1) 
output
}

}
\keyword{ cluster }
