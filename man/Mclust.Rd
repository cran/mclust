\name{Mclust}
\alias{Mclust}
\alias{print.Mclust}

\title{Model-Based Clustering}

\description{
Model-based clustering based on parameterized finite Gaussian mixture models. 
Models are estimated by EM algorithm initialized by hierarchical model-based agglomerative clustering. The optimal model is then selected according to BIC.
}

\usage{
Mclust(data, G = NULL, modelNames = NULL, 
     prior = NULL, 
     control = emControl(), 
     initialization = NULL, 
     warn = mclust.options("warn"), 
     x =  NULL, 
     verbose = interactive(), \dots)
}
\arguments{
  \item{data}{
  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows
  correspond to observations (\eqn{n}) and columns correspond to variables (\eqn{d}). 
  }
  \item{G}{
  An integer vector specifying the numbers of mixture components
  (clusters) for which the BIC is to be calculated. 
  The default is \code{G=1:9}. 
  }
  \item{modelNames}{
  A vector of character strings indicating the models to be fitted 
  in the EM phase of clustering. The default is:
  \itemize{
  \item for univariate data (\eqn{d = 1}): \code{c("E", "V")}
  \item for multivariate data (\eqn{n > d}): all the models available in \code{mclust.options("emModelNames")}
  \item for multivariate data (\eqn{n <= d}): the spherical and diagonal models, i.e. \code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}
  }
  The help file for \code{\link{mclustModelNames}} describes the available models.
   }
  \item{prior}{
  The default assumes no prior, but this argument allows specification of a 
  conjugate prior on the means and variances through the function 
  \code{\link{priorControl}}. \cr
  Note that, as described in \code{\link{defaultPrior}}, in the multivariate 
  case only 10 out of 14 models may be used in conjunction with a prior, i.e.
  those available in \emph{MCLUST} up to version 4.4.
  }
  \item{control}{
  A list of control parameters for EM. The defaults are set by the call
  \code{emControl()}. 
  }
  \item{initialization}{
  A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
  A matrix of merge pairs for hierarchical clustering such as produced
  by function \code{\link{hc}}. \cr
  For multivariate data, the default is to compute a hierarchical 
  agglomerative clustering tree by applying function 
  \code{\link{hc}} with 
  model specified by \code{mclust.options("hcModelName")}, and
  data transformation set by \code{mclust.options("hcUse")}.\cr
  All the input or a subset as indicated by the \code{subset} argument is 
  used for initial clustering.\cr
  The hierarchical clustering results are then used to start the EM
  algorithm from a given partition.\cr
  For univariate data, the default is to use quantiles to start the EM
  algorithm. However, hierarchical clustering could also be used by 
  calling \code{\link{hc}} with model specified as \code{"V"} or \code{"E"}.
  }
  \item{\code{subset}}{
  A logical or numeric vector specifying a subset of the data
  to be used in the initial hierarchical clustering phase.
  By default no subset is used unless the number of observations exceeds 
  the value specified by \code{mclust.options("subset")}. Note that to
  guarantee exact reproducibility of results a seed must be specified 
  (see \code{\link{set.seed}}).}
  \item{\code{noise}}{
  A logical or numeric vector indicating an initial guess as to
  which observations are noise in the data. If numeric the entries
  should correspond to row indexes of the data. If supplied, a noise
  term will be added to the model in the estimation.}
   }
  }
  \item{warn}{
  A logical value indicating whether or not certain warnings 
  (usually related to singularity) should be issued.
  The default is controlled by \code{\link{mclust.options}}.
  }
  \item{x}{
   An object of class \code{'mclustBIC'}. If supplied, BIC values for models 
   that have already been computed and are available in \code{x} are not 
   recomputed. 
   All arguments, with the exception of \code{data}, \code{G} and 
   \code{modelName}, are ignored and their values are set as specified in the
   attributes of \code{x}. Defaults for \code{G} and \code{modelNames}
   are taken from \code{x}.
  }
  \item{verbose}{
  A logical controlling if a text progress bar is displayed during the
  fitting procedure. By default is \code{TRUE} if the session is 
  interactive, and \code{FALSE} otherwise.
  }
  \item{\dots}{
  Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
 An object of class \code{'Mclust'} providing the optimal (according to BIC)
 mixture model estimation.
 
 The details of the output components are as follows:
  
  \item{call}{The matched call} 

  \item{data}{The input data matrix.} 

  \item{modelName}{
  A character string denoting the model at which the optimal BIC occurs.
  }
  \item{n}{
  The number of observations in the data.
  }
  \item{d}{
  The dimension of the data.
  }
  \item{G}{
  The optimal number of mixture components.
  }
  \item{BIC}{
  All BIC values.
  }
  \item{loglik}{
  The log-likelihood corresponding to the optimal BIC.
  }
  \item{df}{
  The number of estimated parameters.
  }
  \item{bic}{
  BIC value of the selected model.
  }
  \item{icl}{
  ICL value of the selected model.
  }
  \item{hypvol}{
  The hypervolume parameter for the noise component if required, otherwise set to \code{NULL} (see \code{\link{hypvol}}).
  }
  \item{parameters}{
   A list with the following components:
   \describe{
    \item{\code{pro}}{
        A vector whose \emph{k}th component is the mixing proportion for
        the \emph{k}th component of the mixture model.
        If missing, equal proportions are assumed.
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
  \item{z}{
  A matrix whose \emph{[i,k]}th entry is the probability that observation
  \emph{i} in the test data belongs to the \emph{k}th class.
  }
  \item{classification}{
  The classification corresponding to \code{z}, i.e. \code{map(z)}.
  }
  \item{uncertainty}{
   The uncertainty associated with the classification.
  }
}
\references{
Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317. 

Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, \emph{Journal of the American Statistical Association}, 97/458, pp. 611-631.

Fraley C., Raftery A. E., Murphy T. B. and Scrucca L. (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation. \emph{Technical Report} No. 597, Department of Statistics, University of Washington.

C. Fraley and A. E. Raftery (2007) Bayesian regularization for normal mixture estimation and model-based clustering. \emph{Journal of Classification}, 24, 155-181.
}
\seealso{
  \code{\link{summary.Mclust}}, 
  \code{\link{plot.Mclust}}, 
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{hc}},
  \code{\link{mclustBIC}}, 
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
mod1 <- Mclust(iris[,1:4])
summary(mod1)

mod2 <- Mclust(iris[,1:4], G = 3)
summary(mod2, parameters = TRUE)

# Using prior
mod3 <- Mclust(iris[,1:4], prior = priorControl())
summary(mod3)

mod4 <- Mclust(iris[,1:4], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(mod4)

# Clustering of faithful data with some artificial noise added 
nNoise <- 100
set.seed(0) # to make it reproducible
Noise <- apply(faithful, 2, function(x) 
              runif(nNoise, min = min(x)-.1, max = max(x)+.1))
data <- rbind(faithful, Noise)
plot(faithful)
points(Noise, pch = 20, cex = 0.5, col = "lightgrey")
set.seed(0)
NoiseInit <- sample(c(TRUE,FALSE), size = nrow(faithful)+nNoise, 
          replace = TRUE, prob = c(3,1)/4)
mod5 <- Mclust(data, initialization = list(noise = NoiseInit))
summary(mod5, parameter = TRUE)
plot(mod5, what = "classification")
}

\keyword{cluster}
