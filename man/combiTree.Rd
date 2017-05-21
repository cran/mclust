\name{combiTree}
\alias{combiTree}

\title{Tree structure obtained from combining mixture components}

\description{The clustering structure 

The method implemented in \code{\link{clustCombi}} can be used for combining Gaussian mixture components for clustering. This provides a hierarchical structure which can be graphically represented as a tree.}

\usage{
combiTree(object, what = c("entropy", "step"), 
                  type = c("triangle", "rectangle"), 
                  edgePar = list(col = "darkgray", lwd = 2), 
                  \dots)
}

\arguments{
  \item{object}{
  An object of class \code{'clustCombi'} resulting from a call to \code{\link{clustCombi}}.
}
  \item{what}{
  A string specifying the quantity used to draw the vertical axis. Possible values are \code{"entropy"} (default), and \code{"step"}.
}
  \item{type}{
  A string specifying the dendrogram's type. Possible values are \code{"triangle"} (default), and \code{"rectangle"}.
}
  \item{edgePar}{
  A list of plotting parameters. See \code{\link[stats]{dendrogram}}.
}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%\details{}
\value{
The function always draw a tree and invisibly returns an object of class \code{'dendrogram'} for fine tuning.
}
%\references{}
\author{L. Scrucca}
%\note{}

\seealso{\code{\link{clustCombi}}}

\examples{
\dontrun{
data(Baudry_etal_2010_JCGS_examples)
output <- clustCombi(data = ex4.1) 
combiTree(output)
combiTree(output, type = "rectangle")
combiTree(output, what = "step")
combiTree(output, what = "step", type = "rectangle")}
}

\keyword{cluster}
\keyword{hplot}
