\name{mclustVariance}
\alias{mclustVariance}
\title{
  Template for variance specification for parameterized 
  Gaussian mixture models.
}
\description{
   Specification of variance parameters for the various types
   of Gaussian mixture models.
}
\details{
    \itemize{
      \item The \code{variance} component in the {parameters} list from the
            output to e.g. \code{me} or\code{mstep} or input to e.g.
           \code{estep}  may contain one or more of the following
            arguments, depending on the model:
      \describe{
        \item{modelName}{
              A character string indicating the model.
        }
        \item{d}{
              The dimension of the data.
        }
        \item{G}{
              The number of components in the mixture model.
        }
	\item{sigmasq}{
	  for the one-dimensional models ("E", "V") and spherical models
	  ("EII", "VII"). This is either a vector whose \emph{k}th
	  component is the variance for the \emph{k}th component in the
	  mixture model ("V" and "VII"), or a scalar giving the common
	  variance for all components in the mixture model ("E" and
	  "EII").
	}
	\item{Sigma}{
	  For the equal variance models "EII", "EEI", and "EEE". 
          A \emph{d} by \emph{d}  matrix giving the common covariance for all  
         components of the  mixture model.
	}
	\item{cholSigma}{
	  For the equal variance model "EEE". 
          A \emph{d} by \emph{d} upper triangular matrix giving the 
          Cholesky factor of the common covariance for all  
          components of the  mixture model.
	}
	\item{sigma}{
	  For all multidimensional mixture models. A
	  \emph{d} by \emph{d} by \emph{G} matrix array whose
	  \code{[,,k]}th entry is the covariance matrix for
	  the \emph{k}th component of the mixture model. 
	}
	\item{cholsigma}{
	  For the unconstrained covaraince mixture model "VVV". 
	  A \emph{d} by \emph{d} by \emph{G} matrix array whose
	  \code{[,,k]}th entry is the upper triangular Cholesky factor
          of the covariance matrix for the \emph{k}th component of the 
          mixture model. 
	}
        \item{scale}{
            For diagonal models "EEI", "EVI", "VEI", "VVI" and
            constant-shape models "EEV" and "VEV".
            Either a \emph{G}-vector giving the scale of the covariance (the
            \emph{d}th root of its determinant) for each component in the
            mixture model, or a single numeric value if the scale is the
            same for each component.
        }
        \item{shape}{
            For diagonal models "EEI", "EVI", "VEI", "VVI" and
            constant-shape models "EEV" and "VEV".
            Either a \emph{G} by \emph{d} matrix in which the \emph{k}th
            column is the shape of the covariance matrix (normalized to have
            determinant 1) for the \emph{k}th component, or a
            \emph{d}-vector giving a common shape for all components.
        }
         \item{orientation}{
              For the constant-shape models "EEV" and "VEV".
              Either a \emph{d} by \emph{d} by \emph{G} array whose
              \code{[,,k]}th entry is the orthonomal matrix of eigenvectors of
              the covariance matrix of the \emph{k}th component, or a \emph{d}
              by \emph{d} orthonormal matrix if the mixture components have a
              common orientation. The \code{orientation} component
              is not needed in spherical and diagonal models, since
              the principal components are parallel to the coordinate axes 
              so that the orientation matrix is the identity.
           }
        }
      }
  In all cases, the value
 \code{-1} is used as a placeholder for unknown nonzero entries. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}.

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\keyword{cluster}
