\name{mclustVariance}
\alias{mclustVariance}

\title{
Template for variance specification for parameterized Gaussian mixture models
}
\description{
   Specification of variance parameters for the various types
   of Gaussian mixture models.
}
\usage{
mclustVariance(modelName, d = NULL, G = 2)
}

\arguments{
  \item{modelName}{A character string specifying the model.}
  \item{d}{A integer specifying the dimension of the data.}
  \item{G}{An integer specifying the number of components in the mixture model.}  
}

\details{The \code{variance} component in the {parameters} list from the
output to e.g. \code{me} or \code{mstep} or input to e.g. \code{estep} may contain one or more of the following arguments, depending on the model:
\describe{
  \item{\code{modelName}}{
        A character string indicating the model.
  }
  \item{\code{d}}{
        The dimension of the data.
  }
  \item{\code{G}}{
        The number of components in the mixture model.
  }
\item{\code{sigmasq}}{
  for the one-dimensional models (\code{"E"}, \code{"V"}) and spherical
  models (\code{"EII"}, \code{"VII"}). This is either a vector whose
  \emph{k}th component is the variance for the \emph{k}th component in 
  the mixture model (\code{"V"} and \code{"VII"}), or a scalar giving 
  the common variance for all components in the mixture model (\code{"E"}
  and \code{"EII"}).
}
\item{\code{Sigma}}{
  For the equal variance models \code{"EII"}, \code{"EEI"}, and
  \code{"EEE"}. 
  A \emph{d} by \emph{d}  matrix giving the common covariance for all  
  components of the  mixture model.
}
\item{\code{cholSigma}}{
  For the equal variance model {"EEE"}. 
  A \emph{d} by \emph{d} upper triangular matrix giving the 
  Cholesky factor of the common covariance for all  
  components of the  mixture model.
}
\item{\code{sigma}}{
  For all multidimensional mixture models. A
  \emph{d} by \emph{d} by \emph{G} matrix array whose
  \code{[,,k]}th entry is the covariance matrix for
  the \emph{k}th component of the mixture model. 
}
\item{\code{cholsigma}}{
  For the unconstrained covariance mixture model \code{"VVV"}. 
  A \emph{d} by \emph{d} by \emph{G} matrix array whose
  \code{[,,k]}th entry is the upper triangular Cholesky factor
  of the covariance matrix for the \emph{k}th component of the 
  mixture model. 
}
\item{\code{scale}}{
    For diagonal models \code{"EEI"}, \code{"EVI"}, \code{"VEI"}, 
    \code{"VVI"} and constant-shape models \code{"EEV"} and \code{"VEV"}.
    Either a \emph{G}-vector giving the scale of the covariance (the
    \emph{d}th root of its determinant) for each component in the
    mixture model, or a single numeric value if the scale is the
    same for each component.
}
\item{\code{shape}}{
    For diagonal models \code{"EEI"}, \code{"EVI"}, \code{"VEI"}, 
    \code{"VVI"} and constant-shape models \code{"EEV"} and \code{"VEV"}.
    Either a \emph{G} by \emph{d} matrix in which the \emph{k}th
    column is the shape of the covariance matrix (normalized to have
    determinant 1) for the \emph{k}th component, or a
    \emph{d}-vector giving a common shape for all components.
}
 \item{\code{orientation}}{
      For the constant-shape models \code{"EEV"} and \code{"VEV"}.
      Either a \emph{d} by \emph{d} by \emph{G} array whose
      \code{[,,k]}th entry is the orthonomal matrix whose
      columns are the eigenvectors of the covariance matrix of
      the \emph{k}th component, or a \emph{d} by \emph{d}
      orthonormal matrix if the mixture components have a
      common orientation. The \code{orientation} component
      is not needed in spherical and diagonal models, since
      the principal components are parallel to the coordinate axes 
      so that the orientation matrix is the identity.
   }
}
In all cases, the value
\code{-1} is used as a placeholder for unknown nonzero entries. 
}

\keyword{cluster}
