\name{em}
\alias{em}

\title{EM algorithm starting with E-step for parameterized Gaussian mixture models}

\description{
  Implements the EM algorithm for parameterized Gaussian mixture models,
  starting with the expectation step.
}
\usage{
em(data, modelName, parameters, prior = NULL, control = emControl(),
   warn = NULL, \dots)
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
     A names list giving the parameters of the model.
     The components are as follows:
     \describe{
         \item{\code{pro}}{
               Mixing proportions for the components of the mixture. 
               If the model includes a Poisson term for noise, there 
               should be one more mixing proportion than the number 
               of Gaussian components.
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
        \item{\code{Vinv}}{
              An estimate of the reciprocal hypervolume of the data region.
              If set to NULL or a negative value, the default is determined by 
              applying function \code{hypvol} to the data.
              Used only when \code{pro} includes an additional
              mixing proportion for a noise component.
        }
     }
  }
 \item{prior}{
      Specification of a conjugate prior on the means and variances.
      The default assumes no prior.
 }
 \item{control}{
       A list of control parameters for EM. The defaults are set by the call 
       \code{emControl()}.
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
A list including the following components: 
  \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{n}{
    The number of observations in the data.
  }
  \item{d}{
    The dimension of the data.
  }
  \item{G}{
    The number of mixture components.
  }  
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{parameters}{
     \describe{
        \item{\code{pro}}{
              A vector whose \emph{k}th component is the mixing proportion for 
              the \emph{k}th component of the mixture model.
              If the model includes a Poisson term for noise, there 
              should be one more mixing proportion than the number 
              of Gaussian components.
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
        \item{\code{Vinv}}{
              The estimate of the reciprocal hypervolume of the data region
              used in the computation when the input indicates the
              addition of a noise component to the model.
        }
      }
  }
  \item{loglik}{
    The log likelihood for the data in the mixture model. 
  }
  \item{control}{
     The list of control parameters for EM used.
  }
  \item{prior}{
     The specification of a conjugate prior on the means and variances used,
     \code{NULL} if no prior is used.
  }        
  \item{Attributes:}{
      \code{"info"} Information on the iteration.\cr
      \code{"WARNING"} An appropriate warning if problems are 
      encountered in the computations.
    }
}
\seealso{
  \code{\link{emE}}, \dots,
  \code{\link{emVVV}},
  \code{\link{estep}},
  \code{\link{me}},
  \code{\link{mstep}},
  \code{\link{mclust.options}},
  \code{\link{do.call}}
}
\examples{
\donttest{
msEst <- mstep(modelName = "EEE", data = iris[,-5], 
               z = unmap(iris[,5]))
names(msEst)

em(modelName = msEst$modelName, data = iris[,-5],
   parameters = msEst$parameters)
}
}
\keyword{cluster}
