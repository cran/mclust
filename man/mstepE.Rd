\name{mstepE}
\alias{mstepE}
\alias{mstepV}
\alias{mstepEII}
\alias{mstepVII}
\alias{mstepEEI}
\alias{mstepVEI}
\alias{mstepEVI}
\alias{mstepVVI}
\alias{mstepEEE}
\alias{mstepEEV}
\alias{mstepVEV}
\alias{mstepVVV}
\alias{mstepEVE}
\alias{mstepEVV}
\alias{mstepVEE}
\alias{mstepVVE}

\title{
  M-step for a parameterized Gaussian mixture model.
}
\description{
  Maximization step in the EM algorithm for a parameterized Gaussian
  mixture model.
}
\usage{
mstepE( data, z, prior = NULL, warn = NULL, \dots)
mstepV( data, z, prior = NULL, warn = NULL, \dots)
mstepEII( data, z, prior = NULL, warn = NULL, \dots)
mstepVII( data, z, prior = NULL, warn = NULL, \dots)
mstepEEI( data, z, prior = NULL, warn = NULL, \dots)
mstepVEI( data, z, prior = NULL, warn = NULL, control = NULL, \dots)
mstepEVI( data, z, prior = NULL, warn = NULL, \dots)
mstepVVI( data, z, prior = NULL, warn = NULL, \dots)
mstepEEE( data, z, prior = NULL, warn = NULL, \dots)
mstepEEV( data, z, prior = NULL, warn = NULL, \dots)
mstepVEV( data, z, prior = NULL, warn = NULL, control = NULL,\dots)
mstepVVV( data, z, prior = NULL, warn = NULL, \dots)
mstepEVE( data, z, prior = NULL, warn = NULL, control = NULL, \dots)
mstepEVV( data, z, prior = NULL, warn = NULL, \dots)
mstepVEE( data, z, prior = NULL, warn = NULL, control = NULL, \dots)
mstepVVE( data, z, prior = NULL, warn = NULL, control = NULL, \dots)
}
\arguments{
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
  \item{control}{
    Values controlling termination for models \code{"VEI"} and \code{"VEV"} 
    that have an iterative M-step. This should be a list with components
    named \emph{itmax} and \emph{tol}. These components can be of length 1 
    or 2; in the latter case, \code{mstep} will use the second value, under 
    the assumption that the first applies to an outer iteration (as in the 
    function \code{me}).
    The default uses the default values from the function \code{emControl},
    which sets no limit on  the number of iterations, and a relative tolerance 
    of \code{sqrt(.Machine$double.eps)} on successive iterates.
   }
  \item{\dots }{
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
       (\code{"VEI"} and \code{"VEV"}), information on the iteration.\cr
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
  \code{\link{mstep}},
  \code{\link{me}},
  \code{\link{estep}},
  \code{\link{mclustVariance}},
  \code{\link{priorControl}},
  \code{\link{emControl}}.
}
\examples{
\dontrun{
mstepVII(data = iris[,-5], z = unmap(iris[,5]))}
}
\keyword{cluster}
