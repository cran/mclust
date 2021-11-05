\name{mclustICL}
\alias{mclustICL}
\alias{print.mclustICL}
\alias{summary.mclustICL}
\alias{print.summary.mclustICL}

\title{ICL Criterion for Model-Based Clustering}

\description{
ICL (Integrated Complete-data Likelihood) for parameterized Gaussian mixture models fitted by EM algorithm initialized by model-based hierarchical clustering.
}

\usage{
mclustICL(data, G = NULL, modelNames = NULL, 
          initialization = list(hcPairs = NULL, 
                                subset = NULL, 
                                noise = NULL), 
          x = NULL, \dots)

\method{summary}{mclustICL}(object, G, modelNames, \dots)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  
  \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the criteria should be calculated. 
    The default is \code{G = 1:9}. 
  }
  
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. The help file for 
    \code{\link{mclustModelNames}} describes the available models.
    The default is:
    \describe{
        \item{\code{c("E", "V")}}{for univariate data}
        \item{\code{mclust.options("emModelNames")}}{for multivariate data (n > d)}
        \item{\code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}}{the spherical and diagonal models for multivariate data (n <= d)}
     }
   }

 \item{initialization}{
    A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. For multivariate data, the default is to compute
    a hierarchical clustering tree by applying function \code{hc} with
    \code{modelName = "VVV"} to the data or a subset as indicated by the
    \code{subset} argument.
    The hierarchical clustering results are to start EM.
    For univariate data, the default is to use quantiles to start EM.
  }
  \item{\code{subset}}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.
  }
  }
  }
  \item{x}{
     An object of class \code{'mclustICL'}. If supplied, \code{mclustICL}
     will use the settings in \code{x} to produce another object of
     class \code{'mclustICL'}, but with \code{G} and \code{modelNames}
     as specified in the arguments. Models that have already been computed
     in \code{x} are not recomputed. All arguments to \code{mclustICL} 
     except \code{data}, \code{G} and \code{modelName} are
     ignored and their values are set as specified in the attributes of
     \code{x}. 
     Defaults for \code{G} and \code{modelNames} are taken from \code{x}.
  }
  \item{\dots}{
    Futher arguments used in the call to \code{\link{Mclust}}. 
    See also \code{\link{mclustBIC}}.
  }

  \item{object}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the criteria should be calculated. 
    The default is \code{G = 1:9}. 
  }

}

\value{
Returns an object of class \code{'mclustICL'} containing the the ICL criterion 
for the specified mixture models and numbers of clusters.

The corresponding \code{print} method shows the matrix of values and the top models according to the ICL criterion. The \code{summary} method shows only the top models.
}

\references{
Biernacki, C., Celeux, G., Govaert, G. (2000). 
Assessing a mixture model for clustering with the integrated completed likelihood.
\emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22 (7), 719-725.

Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317. 
}

\seealso{
  \code{\link{plot.mclustICL}}, 
  \code{\link{Mclust}}, 
  \code{\link{mclustBIC}}, 
  \code{\link{mclustBootstrapLRT}}, 
  \code{\link{bic}},
  \code{\link{icl}}
}

\examples{
data(faithful)
faithful.ICL <- mclustICL(faithful)
faithful.ICL
summary(faithful.ICL)
plot(faithful.ICL)
\donttest{
# compare with
faithful.BIC <- mclustBIC(faithful)
faithful.BIC
plot(faithful.BIC)
}
}

\keyword{cluster}
