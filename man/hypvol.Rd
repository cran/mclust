\name{hypvol}
\alias{hypvol}
\title{
  Aproximate Hypervolume for Multivariate Data
}
\description{
  Computes a simple approximation to the hypervolume of a multivariate
  data set.
}
\usage{
hypvol(data, reciprocal=FALSE)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{reciprocal}{
    A logical variable indicating whether or not the reciprocal
    hypervolume is desired rather than the hypervolume itself. The
    default is to return the hypervolume.  
  }
}
\value{
  Returns the minimum of the hypervolume computed from simple variable bounds
  and that computed from variable bounds of the principal component scores.
  Used for the default hypervolume parameter for the noise 
  component when observations are designated as noise in  \code{Mclust}
  and \code{mclustBIC}.
}
\section{References}{
  A. Dasgupta and A. E. Raftery (1998).
  Detecting features in spatial point processes with clutter via model-based
  clustering. 
 \emph{Journal of the American Statistical Association 93:294-302}. 

  C. Fraley and A.E. Raftery (1998).
   \emph{Computer Journal 41:578-588}.

  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
}
\seealso{
  \code{\link{mclustBIC}} 
}
\examples{
hypvol(iris[,-5])
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
