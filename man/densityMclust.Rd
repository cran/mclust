\name{densityMclust}
\alias{densityMclust}
\title{
  Density Estimation via Model-Based Clustering
}
\description{
  Produces a density estimate for each data point using the optimal
  mixture model from \code{Mclust}.
}
\usage{
densityMclust(x, ...)

}
\arguments{
  \item{x}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{\dots }{
    Additional arguments for the \code{Mclust} function.
  }
}
\value{
  An object of class \code{densityMclust}, which is the density
  of the input $x$ according to the optimal model from \code{Mclust},
  accompanied by the following attributes:
  \item{modelName}{
    A character string denoting the model at which the optimal BIC occurs.
  }
  \item{parameters}{
    The parameters of the \code{Mclust} model.
  }
  \item{range}{
     The range of the input data (used as a default for plotting).
  }
}
\details{
  If \code{x} is univariate, there is an associated \code{plot} method.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  
  C. Fraley and A. E. Raftery (2006, revised in 2010).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics,  University of Washington.
}
\seealso{
  \code{\link{Mclust}}, 
  \code{\link{plot.densityMclust}}, 
}
\examples{

  densityMclust(rivers)
\dontrun{
 plot(densityMclust(rivers))
 plot(densityMclust(rivers), data = rivers, xlim = c(0,1500))
}
}
\keyword{cluster}
% docclass is function
