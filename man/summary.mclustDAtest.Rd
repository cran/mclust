\name{summary.mclustDAtest}
\alias{summary.mclustDAtest}
\title{
  Classification and posterior probability from mclustDAtest.
}
\description{
  Classifications from \code{mclustDAtest} and the corresponding
  posterior probabilities. 
}
\usage{
summary.mclustDAtest(object, pro, \dots)
}
\arguments{
  \item{object}{
    The output of \code{mclustDAtest}.
  }
  \item{pro}{
    Prior probabilities for each class in the training data.
  }
  \item{\dots}{
    Not used. For generic/method consistency.
  }
}
\value{
  A list with the following two components: 
  
  \item{classfication}{
    The classification from \code{mclustDAtest} 
  }
  \item{z}{
    Matrix of posterior probabilities in which the \code{[i,j]}th entry
    is the probability of observation \emph{i} belonging to class
    \emph{j}.
  }
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}.
  See \url{http://www.stat.washington.edu/mclust}. 
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington.
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{mclustDAtest}}
}
\examples{
set.seed(0)
n <- 100 ## create artificial data

x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))
\dontrun{
par(pty = "s")
mclust2Dplot(x, classification = xclass, type="classification", ask=FALSE)
}

odd <- seq(1, 2*n, 2)
train <- mclustDAtrain(x[odd, ], labels = xclass[odd]) ## training step
summary(train)

even <- seq(1, 2*n, 2)
test <- mclustDAtest(x[even, ], train) ## compute model densities
testSummary <- summary(test) ## classify training set

names(testSummary)
testSummary$class
testSummary$z
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
