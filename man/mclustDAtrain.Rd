\name{mclustDAtrain}
\alias{mclustDAtrain}
\alias{print.mclustDAtrain}
\title{
  MclustDA Training
}
\description{
  Training phase for MclustDA discriminant analysis.
}
\usage{
mclustDAtrain(data, labels, G=NULL, modelNames=NULL, prior=NULL, 
              control=emControl(), initialization=NULL, warn=FALSE, 
              verbose=TRUE, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{labels}{
    A numeric or character vector assigning a class label to each observation.
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
    A logical value indicating whether or not to print the models and
    numbers of components for each class.
    Default: \code{verbose=TRUE}.
  }
  \item{\dots }{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A list in which each element gives the parameters and other summary
  information for the model best fitting each class according to BIC.
  Attributes are the input parameters other than \code{data}, \code{labels}
  and \code{verbose}.
}
\details{
  Except for \code{labels} and \code{verbose}, the arguments are the
  same as those for \code{mclustBIC}.
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
  \code{\link{summary.mclustDAtrain}}, 
  \code{\link{mclustDAtest}}, 
  \code{\link{mclustBIC}} 
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
