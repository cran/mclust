\name{mclustDA}
\alias{mclustDA}
\alias{print.mclustDA}
\title{
  MclustDA discriminant analysis.
}
\description{
  MclustDA training and testing.
}
\usage{
mclustDA(trainingData, labels, testData, G=1:6, verbose = FALSE)
}
\arguments{
  \item{trainingData}{
    
    A numeric vector, matrix, or data frame of training observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{labels}{
    A numeric or character vector assigning a class label to each 
    training observation.
  }
  \item{testData}{
    
    A numeric vector, matrix, or data frame of training observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{G}{
    An integer vector specifying the numbers of mixture components (clusters)
    to be considered for each class.
    Default: \code{1:6}.
  }
  \item{verbose}{
    A logical variable telling whether or not to print an indication that the
    function is in the training phase, which may take some time to complete.
  }
}
\value{
  A list with the following components: 
  \item{testClassification}{ 
    \code{mclustDA} classification of the test data.
  }
  \item{trainingClassification}{ 
    \code{mclustDA} classification of the training data.
  }
  \item{errorRate}{ 
    The error rate for the \code{mclustDA} classification 
    of the training data compared to the known labels,
  }
  \item{summary}{ 
    Gives the best model and number of clusters for each training class.
  }
  \item{models}{ 
    The mixture models used to fit the known classes.
  }
  \item{postProb}{ 
    A  matrix whose \emph{[i,k]}th entry is the probability that observation
    \emph{i} in the test data belongs to the \emph{k}th class.
  }
}
\section{Details}{
  The following models are compared in \code{Mclust}: \cr\cr
  "E" for spherical, equal variance (one-dimensional) \cr
  "V" for spherical, variable variance (one-dimensional) \cr\cr
  "EII": spherical, equal volume \cr
  "VII": spherical, unequal volume \cr
  "EEI": diagonal, equal volume, equal shape \cr
  "VVI": diagonal, varying volume, varying shape \cr
  "EEE": ellipsoidal, equal volume, shape, and orientation \cr
  "VVV": ellipsoidal, varying volume, shape, and orientation 

  \code{mclustDA} is a simplified function combining
  \code{mclustDAtrain} and \code{mclustDAtest} and their summaries.
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
  \code{\link{plot.mclustDA}}, 
  \code{\link{mclustDAtrain}}, 
  \code{\link{mclustDAtest}}
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
odd <- seq(from = 1, to = 2*n, by = 2)
even <- odd + 1
testMclustDA <- mclustDA(trainingData = x[odd, ], labels = xclass[odd], 
                         testData = x[even,])

clEven <- testMclustDA$testClassification ## classify training set
compClass(clEven,xclass[even])$error

\dontrun{plot(testMclustDA, trainingData = x[odd, ], labels = xclass[odd], 
              testData = x[even,])}
}
\keyword{cluster}
% docclass is function
