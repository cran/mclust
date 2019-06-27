\name{mstep}
\alias{mstep}

\title{M-step for parameterized Gaussian mixture models}

\description{
  Maximization step in the EM algorithm for parameterized Gaussian
  mixture models.
}
\usage{
mstep(modelName, data, z, prior = NULL, warn = NULL, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{z}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the ith observation belonging to
    the \emph{k}th component of the mixture.  
    In analyses involving noise, this should not include the
    conditional probabilities for the noise component. 
  }
 \item{prior}{                                                      
     Specification of a conjugate prior on the means and variances.
     The default assumes no prior.
   }
 \item{warn}{
    A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when the
    estimation fails. The default is given by \code{mclust.options("warn")}.
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
      }
  }
  \item{Attributes:}{
    \code{"info"} For those models with iterative M-steps
       (\code{"VEI"} and \code{"VEV"}), information  on the iteration.\cr
    \code{"WARNING"} An appropriate warning if problems are
    encountered in the computations.
  }
}
\note{
   This function computes the M-step only for MVN mixtures, so in 
   analyses involving noise, the conditional probabilities input should 
   exclude those for the noise component. \cr

   In contrast to \code{me} for the EM algorithm, computations in \code{mstep}
   are carried out unless failure due to overflow would occur. To impose
   stricter tolerances on a single \code{mstep}, use \code{me} with the
  \emph{itmax} component of the \code{control} argument set to 1.
}

\seealso{
  \code{\link{mstepE}}, \dots,
  \code{\link{mstepVVV}},
  \code{\link{emControl}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclust.options}}.
}
\examples{
\dontrun{
mstep(modelName = "VII", data = iris[,-5], z = unmap(iris[,5]))}
}
\keyword{cluster}
