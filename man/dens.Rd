\name{dens}
\alias{dens}
\title{
  Density for Parameterized MVN Mixtures
}
\description{
  Computes densities of observations in parameterized MVN mixtures.
}
\usage{
dens(modelName, data, logarithm = FALSE, parameters, warn=NULL, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{logarithm}{
     A logical value indicating whether or not the logarithm of the component 
      densities should be returned. The default is to return the component 
      densities, obtained from the log component densities by exponentiation.
   }
 \item{parameters}{
     The parameters of the model:
     \describe{ 
       \item{\code{pro}}{
           The vector of mixing proportions for the components of the mixture. 
        }
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th 
              component of the mixture model. 
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}} 
              for details.  
            }
 
     }
  }
  \item{warn}{
    A logical value indicating whether or not a warning should be issued
    when computations fail. The default is \code{warn=FALSE}.
  }
  \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A numeric vector whose \emph{i}th component is the density of the
  \emph{ith} observation in \code{data} in the MVN mixture specified 
  by \code{parameters}.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{cdens}},
  \code{\link{mclust.options}},
  \code{\link{do.call}}
}
\examples{
faithfulModel <- Mclust(faithful)
Dens <- dens(modelName = faithfulModel$modelName, data = faithful,
             parameters = faithfulModel$parameters)
Dens

\dontrun{
## alternative call
do.call("dens", faithfulModel)
}
}
\keyword{cluster}

