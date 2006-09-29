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
mclustDA(train, test, pro=NULL, G=NULL, modelNames=NULL, prior=NULL, 
         control=emControl(), initialization=NULL, 
         warn=FALSE, verbose=FALSE, \dots)
}
\arguments{
  \item{train}{
    A list with two named components: \code{data} giving the data
    and \code{labels} giving the class labels for the observations in
    the data.
  }
  \item{test}{
    A list with two named components: \code{data} giving the data
    and \code{labels} giving the class labels for the observations in
    the data. The labels are used only to compute the error rate
    in the \code{print} method and  can be set to \code{NULL} if unknown.
    The default is to test the training data.
  }
 \item{pro}{
    Optional prior probabilities for each class in the training data.
  }
 \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the BIC is to be calculated. 
    The default is \code{G=1:9}. 
  }
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. The help file for
    \code{mclustModelNames} describes the available models.
    The default is \code{c("E", "V")} for univariate data and
    \code{mclustOptions()\$emModelNames} for multivariate data.
   }
  \item{prior}{
    The default assumes no prior, but this argument allows specification of a 
    conjugate prior on the means and variances through the function 
    \code{priorControl}.
    }
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
  \item{initialization}{
   A list containing zero or more of the following components:
  \itemize{
  \item{hcPairs}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. The default is to compute a hierarchical
    clustering tree by applying function \code{hc} with
    \code{modelName = "E"} to univariate data and
    \code{modelName = "VVV"} to multivariate data or a
    subset as indicated by the \code{subset} argument. 
    The hierarchical clustering results are used as starting values for EM.  
  }
  \item{subset}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.
  }
  }
  }
  \item{warn}{
     A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when
     estimation fails. The default is to suppress these warnings.
  }
  \item{verbose}{
    A logical variable telling whether or not to print an indication that the
    function is in the training phase, which may take some time to complete.
  }
  \item{\dots }{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A list with the following components: 
  \item{test}{ 
    A list with the following components:
     \describe{
        \item{classification}{
              The classification of the test data for this instance of 
              \code{mclustDA}.
        }
        \item{uncertainty}{
            The uncertainty of the classification (0 least certain,
            1 most certain).
        }
        \item{labels}{
              The test labels (if any) from the input.
        }
      }
  }
  \item{training}{ 
    A list with the following components:
     \describe{
        \item{classification}{
              The classification of the training data for this instance of 
              \code{mclustDA}.
        }
        \item{z}{
            A  matrix whose \emph{[i,k]}th entry is the probability that 
            observation \emph{i} in the training data belongs to the 
            \emph{k}th class.
        }
        \item{labels}{
              The training labels from the input.
        }
      }
  }
  \item{summary}{ 
    A data frame summarizing the \code{mclustDA} results including
    the mixture models and numbers of components for the training classes.
  }
}
\details{
  \code{mclustDA} combines functions \code{mclustDAtrain} and 
  \code{mclustDAtest} and their summaries. This is suitable when 
   all test data are available in advance, so that the training
   model is only used once.
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
  \code{\link{plot.mclustDA}}, 
  \code{\link{mclustDAtrain}}, 
  \code{\link{mclustDAtest}},
  \code{\link{classError}}
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
