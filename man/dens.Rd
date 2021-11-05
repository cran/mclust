\name{dens}
\alias{dens}
\title{
  Density for Parameterized MVN Mixtures
}
\description{
  Computes densities of observations in parameterized MVN mixtures.
}
\usage{
dens(data, modelName, parameters, logarithm = FALSE, warn=NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
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
  \item{logarithm}{
     A logical value indicating whether or not the logarithm of the component 
      densities should be returned. The default is to return the component 
      densities, obtained from the log component densities by exponentiation.
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

\seealso{
  \code{\link{cdens}},
  \code{\link{mclust.options}},
  \code{\link{do.call}}
}
\examples{
\donttest{
faithfulModel <- Mclust(faithful)
Dens <- dens(modelName = faithfulModel$modelName, data = faithful,
             parameters = faithfulModel$parameters)
Dens

## alternative call
do.call("dens", faithfulModel)}
}
\keyword{cluster}

