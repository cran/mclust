\name{cdens}
\alias{cdens}
\title{
  Component Density for Parameterized MVN Mixture Models
}
\description{
  Computes component densities for observations in MVN mixture models
  parameterized by eigenvalue decomposition.
}
\usage{
cdens(modelName, data, logarithm = FALSE, parameters, warn = NULL, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{mclustModelNames} describes the available models.
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
        \item{mean}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th 
              component of the mixture model. 
        }
        \item{variance}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{mclustVariance} 
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
  A numeric matrix whose \code{[i,k]}th entry is the 
  density or log density of observation \emph{i} in component \emph{k}. 
  The densities are not scaled by mixing proportions.
}
\note{
  When one or more component densities are very large in magnitude,
  it may be possible to compute the logarithm of the component
  densities but not the component densities themselves due to overflow.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics, University of Washington.
}

\seealso{
  \code{\link{cdensE}}, \ldots,
  \code{\link{cdensVVV}},
  \code{\link{dens}},
  \code{\link{estep}},
  \code{\link{mclustModelNames}},
  \code{\link{mclustVariance}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{
z2 <- unmap(hclass(hcVVV(faithful),2)) # initial value for 2 class case

model <- me( modelName="EEE", data=faithful, z=z2)
cdens(modelName="EEE", data=faithful, logarithm = TRUE, 
      parameters = model$parameters)[1:5,]

odd <- seq(1, nrow(cross), by = 2)
oddBIC <- mclustBIC(cross[odd,-1]) 
oddModel <- mclustModel(cross[odd,-1], oddBIC) ## best parameter estimates
names(oddModel)

even <- odd + 1
densities <- cdens(modelName = oddModel$modelName, data = cross[even,-1], 
                   parameters = oddModel$parameters)
cbind(class = cross[even,1], densities)[1:5,]
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
