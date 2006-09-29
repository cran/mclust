\name{summary.mclustModel}
\alias{summary.mclustModel}
\title{
  Summary Function for MCLUST Models
}
\description{
  Classification and uncertainty for a mixture models 
  as output by \code{mclustModel}.
}
\usage{
\method{summary}{mclustModel}(object, \dots)
}
\arguments{
  \item{object}{
    An \code{"mclustModel"} object.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
   A data frame giving the classification and uncertainty 
   corresponding to the model.
 }
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}.

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{mclustModel}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
irisModel <- mclustModel(iris[,-5], irisBIC)
summary(irisModel)
}
\keyword{cluster}
% docclass is function
