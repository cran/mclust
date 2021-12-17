\name{crimcoords}
\alias{crimcoords}
\alias{print.crimcoords}
\alias{plot.crimcoords}

\title{Discriminant coordinates data projection}

\description{
Compute the discriminant coordinates or crimcoords obtained by projecting the observed data from multiple groups onto the discriminant subspace.
The optimal projection subspace is given by the linear transformation of the original variables that maximizes the ratio of the between-groups covariance (which represents groups separation) to the pooled within-group covariance (which represents within-group dispersion).}

\usage{
crimcoords(data, classification, 
           numdir = NULL, 
           unbiased = FALSE, 
           plot = TRUE, \dots)
           
\method{plot}{crimcoords}(x, \dots)
}

\arguments{
  
  \item{data}{
  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows
  correspond to observations and columns correspond to variables. 
  }
  
  \item{classification}{
	A vector (numerical, character string, or factor) giving the 
	groups classification (either the known class labels or the estimated 
	clusters) for the observed data.}

  \item{numdir}{
	An integer value specifying the number of directions of the 
	discriminant subspace to return. If not provided, the maximal number of
	directions are returned (which is given by the number of non-null 
	eigenvalues, the minimum among the number of variables and the number
	of groups minus one). 
	However, since the effectiveness of the discriminant coordinates in 
	highlighting the separation of groups is decreasing, it might be useful
	to provide a smaller value, say 2 or 3.}

  \item{unbiased}{
	A logical specifying if unbiased estimates should be used for the 
	between-groups and within-groups covariances. By default
	\code{unbiased = FALSE} so MLE estimates are used.
	Note that the use of unbiased or MLE estimates only changes the 
	eigenvalues and eigenvectors of the generalized eigendecomposition by 
	a constant of proportionality, so the discriminant coordinates or 
	crimcoords are essentially the same.}

  \item{plot}{
	A logical specifying if a graph of discriminant coordinates should be 
	produced (default) or not.}
	
	\item{x}{
	An object of class \code{crimcoords} as returned by \code{crimcoords()} function.}
	
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
A list of class \code{crimcoords} with the following components:

\item{means}{A matrix of within-groups means.}
\item{B}{The between-groups covariance matrix.}
\item{W}{The pooled within-groups covariance matrix.}
\item{evalues}{A vector of eigenvalues.}
\item{basis}{A matrix of eigenvectors specifying the basis of the 
discriminant subspace.}
\item{projection}{A matrix of projected data points onto the discriminant
subspace.}
\item{classification}{A vector giving the groups classification.}

If \code{plot = TRUE} (default) a graph of data points projected onto 
the discriminant coordinate or crimcoords subspace is also produced.
}

\references{
Gnanadesikan, R. (1977) \emph{Methods for Statistical Data Analysis of Multivariate Observations}. John Wiley 1& Sons, Sec. 4.2.

Flury, B. (1997) \emph{A First Course in Multivariate Statistics}. Springer, Sec. 7.3.
}

\author{
Luca Scrucca \email{luca.scrucca@unipg.it}
}

%\note{}

\seealso{\code{\link{MclustDR}}, \code{\link{clPairs}}.}

\examples{
# discriminant coordinates for the iris data using known classes 
data("iris")
CRIMCOORDS = crimcoords(iris[,-5], iris$Species)
CRIMCOORDS

# banknote data
data("banknote")

# discriminant coordinate on known classes 
CRIMCOORDS = crimcoords(banknote[,-1], banknote$Status)
CRIMCOORDS

#  discriminant coordinates on estimated clusters
mod = Mclust(banknote[,-1])
CRIMCOORDS = crimcoords(banknote[,-1], mod$classification, plot = FALSE)
plot(CRIMCOORDS$projection, type = "n")
text(CRIMCOORDS$projection, cex = 0.8,
     labels = strtrim(banknote$Status, 2), 
     col = mclust.options("classPlotColors")[1:mod$G][mod$classification])
}

\keyword{multivariate}
