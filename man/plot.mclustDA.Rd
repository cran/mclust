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
plot.mclustDA(x, trainingData, labels, testData, dimens=c(1,2),
              scale = FALSE, identify=FALSE, \dots)
}
\arguments{
  \item{x}{
    The object produced by applying \code{mclustDA} with
    \code{trainingData} and classification \code{labels} to
    \code{testData}. 
  }
  \item{trainingData}{
    The numeric vector, matrix, or data frame of training observations
    used to obtain \code{x}.  
  }
  \item{labels}{
    The numeric or character vector assigning a class label to each
    training observation. 
  }
  \item{testData}{
    A numeric vector, matrix, or data frame of training
    observations. Categorical variables are not allowed. If a matrix or
    data frame, rows correspond to observations and columns correspond
    to variables. 
  }
  \item{dimens}{
    An integer vector of length two specifying the dimensions for coordinate
    projections if the data is more than two-dimensional.
    The default is \code{c(1,2)} (the first two dimesions).
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and thus preserve
    the shape of the distribution. Default: \code{scale=FALSE}  
  }
  \item{identify}{
    A logical variable indicating whether or not to print a title identifying
    the plot. Default: \code{identify=FALSE} 
  }
  \item{\dots}{
    Further arguments to the lower level plotting functions.
  }
}
\value{
  Plots selected via a menu including the following options: 
  training and test data, known training data classification,
  mclustDA test data classification, training errors.
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
  \code{\link{mclustDA}}
}
\examples{
n <- 250 ## create artificial data
set.seed(0)
x <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
xclass <- c(rep(1,n),rep(2,n))
\dontrun{
mclust2Dplot(x, classification = xclass, type="classification", ask=FALSE)
}
odd <- seq(from = 1, to = 2*n, by = 2)
even <- odd + 1
testMclustDA <- mclustDA(trainingData = x[odd, ], labels = xclass[odd], 
                         testData = x[even,])

clEven <- testMclustDA$testClassification ## classify training set
compClass(clEven,xclass[even])$error

\dontrun{
plot(testMclustDA, trainingData = x[odd, ], labels = xclass[odd], 
testData = x[even,])
}
} 
\keyword{cluster}
% docclass is function
