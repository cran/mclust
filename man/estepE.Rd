\name{estepE}
\alias{estepE}
\alias{estepV}
\alias{estepEII}
\alias{estepVII}
\alias{estepEEI}
\alias{estepVEI}
\alias{estepEVI}
\alias{estepVVI}
\alias{estepEEE}
\alias{estepEEV}
\alias{estepVEV}
\alias{estepVVV}
\alias{estepEVE}
\alias{estepEVV}
\alias{estepVEE}
\alias{estepVVE}

\title{
  E-step in the EM algorithm for a parameterized Gaussian mixture model.
}
\description{
  Implements the expectation step in the EM algorithm for a 
  parameterized Gaussian mixture model.
}
\usage{
estepE(data, parameters, warn = NULL, \dots)
estepV(data, parameters, warn = NULL, \dots)
estepEII(data, parameters, warn = NULL, \dots)
estepVII(data, parameters, warn = NULL, \dots)
estepEEI(data, parameters, warn = NULL, \dots)
estepVEI(data, parameters, warn = NULL, \dots)
estepEVI(data, parameters, warn = NULL, \dots)
estepVVI(data, parameters, warn = NULL, \dots)
estepEEE(data, parameters, warn = NULL, \dots)
estepEEV(data, parameters, warn = NULL, \dots)
estepVEV(data, parameters, warn = NULL, \dots)
estepVVV(data, parameters, warn = NULL, \dots)
estepEVE(data, parameters, warn = NULL, \dots)
estepEVV(data, parameters, warn = NULL, \dots)
estepVEE(data, parameters, warn = NULL, \dots)
estepVVE(data, parameters, warn = NULL, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{parameters}{
    The parameters of the model:
    %\itemize{
    %\item An argument describing the variance (depends on the model):
      \describe{
         \item{\code{pro}}{
               Mixing proportions for the components of the mixture.
               If the model includes a Poisson term for noise, there
               should be one more mixing proportion than the number
               of Gaussian components.
        }
        \item{mu}{
              The mean for each component. If there is more than one component,
              this is a matrix whose columns are the means of the  components.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}} 
              for details.
        }
        \item{\code{Vinv}}{
              An estimate of the reciprocal hypervolume of the data region.
              If not supplied or set to a negative value, the default is
              determined by applying function \code{hypvol} to the data.
              Used only when \code{pro} includes an additional
              mixing proportion for a noise component.
        }
      }
    %}
  }
  \item{warn}{
    A logical value indicating whether or certain warnings should be issued.
    The default is given by \code{mclust.options("warn")}.
  }
  \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A list including the following components: 
  \item{modelName}{
    Character string identifying the model.
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the \emph{i}th observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{parameters}{
    The input parameters.
  }
  \item{loglik}{
    The logliklihood for the data in the mixture model. 
  }
  \item{Attribute}{
      \code{"WARNING"}: An appropriate warning if problems are
      encountered in the computations.
  }
}
\seealso{
  \code{\link{estep}},
  \code{\link{em}},
  \code{\link{mstep}},
  \code{\link{do.call}},
  \code{\link{mclustVariance}},
  \code{\link{mclust.options}}.
}
\examples{
\donttest{
msEst <- mstepEII(data = iris[,-5], z = unmap(iris[,5]))
names(msEst)

estepEII(data = iris[,-5], parameters = msEst$parameters)}
}
\keyword{cluster}
