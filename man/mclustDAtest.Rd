\name{mclustDAtest}
\alias{mclustDAtest}
\title{
  MclustDA Testing
}
\description{
  Testing phase for MclustDA discriminant analysis.
}
\usage{
mclustDAtest(data, models)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations
    to be classified.
  }
  \item{models}{
    A list of MCLUST-style models including parameters, usually the
    result of applying \code{mclustDAtrain} to some training data. 
  }
}
\value{
  A matrix in which the \code{[i,j]}th entry is the 
  density for that test observation \emph{i} in the model for class \emph{j}.
}
\details{
  Apply \code{summary} to the output to obtain the classification 
  of the test data.
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
  \code{\link{summary.mclustDAtest}},
  \code{\link{classError}},
  \code{\link{mclustDAtrain}}
}
\examples{
odd <- seq(1, nrow(cross), by = 2)
train <- mclustDAtrain(cross[odd,-1], labels = cross[odd,1]) ## training step
summary(train)

even <- odd + 1
test <- mclustDAtest(cross[even,-1], train) ## compute model densities
clEven <- summary(test)$class ## classify training set
classError(clEven,cross[even,1])
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
