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
  \code{\link{summary.mclustDAtest}},
  \code{\link{mclustDAtrain}}
}
\examples{
n <- 250 ## create artificial data
set.seed(0)
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

even <- odd + 1
test <- mclustDAtest(x[even, ], train) ## compute model densities
summary(test)$class ## classify training set
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
