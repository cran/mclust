\name{plot.mclustDA}
\alias{plot.mclustDA}
\title{
  Plotting method for MclustDA discriminant analysis.
}
\description{
  Plots training and test data, known training data classification,
  mclustDA test data classification, and/or training errors.
}
\usage{
plot.mclustDA(x, trainData, testData, \dots)
}
\arguments{
  \item{x}{
    The object produced by applying \code{mclustDA} with
    \code{trainingData} and classification \code{labels} to
    \code{testData}. 
  }
  \item{trainData}{
    The numeric vector, matrix, or data frame of training observations
    used to obtain \code{x}.  
  }
  \item{testData}{
    A numeric vector, matrix, or data frame of training
    observations. Categorical variables are not allowed. If a matrix or
    data frame, rows correspond to observations and columns correspond
    to variables. 
  }
  \item{\dots}{
    Further arguments to the lower level plotting functions.
  }
}
\value{
  Plots of the following: 
  training and test data, known training data classification,
  mclustDA test data classification, and (if test labels were
  supplied to \code{mclustDA} when \code{x} was created) test errors.
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
  \code{\link{mclustDA}}
}
\examples{
n <- 250 ## create artificial data
set.seed(1)
triModal <- c(rnorm(n,-5), rnorm(n,0), rnorm(n,5))
triClass <- c(rep(1,n), rep(2,n), rep(3,n))

odd <- seq(from = 1, to = length(triModal), by = 2)
even <- odd + 1
triMclustDA <- mclustDA(train=list(data=triModal[odd],labels=triClass[odd]),
                   test= list(data=triModal[even],labels=triClass[even]),
                       verbose = TRUE)

names(triMclustDA)
\dontrun{
  plot(triMclustDA, trainData = triModal[odd], testData = triModal[even])
}

odd <- seq(from = 1, to = nrow(cross), by = 2)
even <- odd + 1
crossMclustDA <- mclustDA( train=list(data=cross[odd,-1],
                                      labels=cross[odd,1]),
                       test= list(data=cross[even,-1],labels=cross[even,1]),
                       verbose = TRUE)

\dontrun{
  plot(crossMclustDA, trainData = cross[odd,-1], testData = cross[even,-1])
}

odd <- seq(from = 1, to = nrow(iris), by = 2)
even <- odd + 1
irisMclustDA <- mclustDA(train=list(data=iris[odd,-5],labels=iris[odd,5]),
                       test= list(data=iris[even,-5],labels=iris[even,5]),
                       verbose = TRUE)

\dontrun{
  plot(irisMclustDA, trainData = iris[odd,-5], testData = iris[even,-5])
}
} 
\keyword{cluster}
% docclass is function
