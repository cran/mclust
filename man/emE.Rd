\name{emE}
\alias{emE}
\alias{emV}
\alias{emX}
\alias{emEII}
\alias{emVII}
\alias{emEEI}
\alias{emVEI}
\alias{emEVI}
\alias{emVVI}
\alias{emEEE}
\alias{emEEV}
\alias{emVEV}
\alias{emVVV}
\alias{emEVV}
\alias{emEVE}
\alias{emVEE}
\alias{emVVE}
\alias{emXII}
\alias{emXXI}
\alias{emXXX}

\title{EM algorithm starting with E-step for a parameterized Gaussian mixture model}

\description{
  Implements the EM algorithm for a parameterized Gaussian mixture model,
  starting with the expectation step.
}
\usage{
emE(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emV(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emX(data, prior = NULL, warn = NULL, \dots)
emEII(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVII(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEEI(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVEI(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEVI(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVVI(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEEE(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVEE(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEVE(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVVE(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEEV(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVEV(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emEVV(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emVVV(data, parameters, prior = NULL, control = emControl(), warn = NULL, \dots)
emXII(data, prior = NULL, warn = NULL, \dots)
emXXI(data, prior = NULL, warn = NULL, \dots)
emXXX(data, prior = NULL, warn = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables.
  }
 \item{parameters}{
     The parameters of the model:
     \describe{
        \item{\code{pro}}{
              Mixing proportions for the components of the mixture. 
              There should one more mixing proportion than the number of 
              Gaussian components if the mixture model includes 
              a Poisson noise term.
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
             The default is determined by applying function  \code{hypvol} 
             to the data. Used only when \code{pro} includes an additional
              mixing proportion for a noise component.
        }
     }
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
  \item{warn}{
    A logical value indicating whether or not a warning should be issued
    whenever a singularity is encountered.
    The default is given in \code{mclust.options("warn")}.
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
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{parameters}{
     \describe{
        \item{\code{pro}}{
              A vector whose \emph{k}th component is the mixing proportion
              for the \emph{k}th component of the mixture model.
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
 \item{Attributes:}{
      \code{"info"} Information on the iteration.\cr
      \code{"WARNING"} An appropriate warning if problems are 
      encountered in the computations.
    }
}
\seealso{
  \code{\link{me}},
  \code{\link{mstep}},
  \code{\link{mclustVariance}},
  \code{\link{mclust.options}}.
}
\examples{
\dontrun{
msEst <- mstepEEE(data = iris[,-5], z = unmap(iris[,5]))
names(msEst)

emEEE(data = iris[,-5], parameters = msEst$parameters)}
}
\keyword{cluster}
