\name{estep}
\alias{estep}
\title{
  E-step for parameterized Gaussian mixture models.
}
\description{
  Implements the expectation step of EM algorithm for parameterized Gaussian
  mixture models. 
}
\usage{
  estep(data, modelName, parameters, warn = NULL, \dots)
}
\arguments{
 \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
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
              If set to NULL or a negative value, the default is determined
              by applying function \code{hypvol} to the data.
              Used only when \code{pro} includes an additional
              mixing proportion for a noise component.
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
  A list including the following components: 
  \item{modelName}{
    A character string identifying the model (same as the input argument).
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the conditional probability
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  }
  \item{parameters}{
    The input parameters.
  }
  \item{loglik}{
    The log-likelihood for the data in the mixture model. 
  }
  \item{Attributes}{
      \code{"WARNING"}: an appropriate warning if problems are
      encountered in the computations.
    }
}
\seealso{
  \code{\link{estepE}}, \dots,
  \code{\link{estepVVV}},
  \code{\link{em}},
  \code{\link{mstep}},
  \code{\link{mclust.options}}
  \code{\link{mclustVariance}}
}
\examples{
\donttest{
msEst <- mstep(modelName = "VVV", data = iris[,-5], z = unmap(iris[,5]))
names(msEst)

estep(modelName = msEst$modelName, data = iris[,-5],
      parameters = msEst$parameters)}
}
\keyword{cluster}
