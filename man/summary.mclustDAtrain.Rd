\name{summary.mclustDAtrain}
\alias{summary.mclustDAtrain}
\title{
  Models and classifications from mclustDAtrain
}
\description{
  The models selected in \code{mclustDAtrain}
  and the corresponding classfications.
}
\usage{
summary.mclustDAtrain(object, \dots)
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
  \code{\link{mclustDAtrain}}
}
\examples{
set.seed(0)
n <- 100 ## create artificial data

x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))
\dontrun{
par(pty = "s")
mclust2Dplot(x, classification = xclass, type="classification", ask=FALSE)}

odd <- seq(1, 2*n, 2)
train <- mclustDAtrain(x[odd, ], labels = xclass[odd]) ## training step
summary(train)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
