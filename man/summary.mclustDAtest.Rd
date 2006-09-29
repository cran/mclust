\name{summary.mclustDAtest}
\alias{summary.mclustDAtest}
\title{
  Classification and posterior probability from mclustDAtest.
}
\description{
  Extract classifications and the corresponding posterior probabilities 
  from \code{mclustDAtest}.
  
}
\usage{
\method{summary}{mclustDAtest}(object, pro=NULL, \dots)
}
\arguments{
  \item{object}{
    The output of \code{mclustDAtest}.
  }
  \item{pro}{
    Optional prior probabilities for each class in the training data.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
  A list with the following two components: 
  
  \item{classfication}{
    The classification from \code{mclustDAtest}.
  }
  \item{z}{
    Matrix of posterior probabilities in which the \code{[i,j]}th entry
    is the probability of observation \emph{i} belonging to class
    \emph{j}.
  }
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
  \code{\link{classError}},
  \code{\link{mclustDAtest}}
}
\examples{
odd <- seq(1, nrow(cross), by = 2)
train <- mclustDAtrain(cross[odd,-1], labels = cross[odd,1]) ## training step
summary(train)

even <- odd + 1
test <- mclustDAtest(cross[even,-1], train) ## compute model densities
testSummary <- summary(test) 
names(testSummary)
classError(testSummary$classification,cross[even,1])
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
