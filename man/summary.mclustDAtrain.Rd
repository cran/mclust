\name{summary.mclustDAtrain}
\alias{summary.mclustDAtrain}
\title{
  Models and classifications from mclustDAtrain
}
\description{
  Extracts the models selected in \code{mclustDAtrain}
  and the corresponding classfications.
}
\usage{
\method{summary}{mclustDAtrain}(object, \dots)
}
\arguments{
  \item{object}{
    The output of \code{mclustDAtrain}.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
  A list identifying the model selected by 
  \code{mclustDAtrain} for each
  class of training data and the corresponding classification.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{mclustDAtrain}}
}
\examples{
odd <- seq(1, nrow(cross), by = 2)
train <- mclustDAtrain(cross[odd,-1], labels = cross[odd,1]) 
summary(train)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
