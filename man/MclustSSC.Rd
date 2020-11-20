\name{MclustSSC}
\alias{MclustSSC}
\alias{print.MclustSSC}

\title{MclustSSC semi-supervised classification}

\description{
Semi-Supervised classification based on Gaussian finite mixture modeling.
}

\usage{
MclustSSC(data, class, 
          G = NULL, modelNames = NULL, 
          prior = NULL, control = emControl(), 
          warn = mclust.options("warn"), 
          verbose = interactive(),
          \dots)
}

\arguments{
  \item{data}{
    A data frame or matrix giving the training data.
  }
  \item{class}{
    A vector giving the known class labels (either a numerical value or 
    a character string) for the observations in the training data. 
    Observations with unknown class are encoded as \code{NA}.
  }
 \item{G}{
    An integer value specifying the numbers of mixture components or classes. 
    By default is set equal to the number of known classes.
    See the examples below.
  }
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    by EM (see the description in \code{\link{mclustModelNames}}).
    See the examples below.
  }
  \item{prior}{
    The default assumes no prior, but this argument allows specification of a 
    conjugate prior on the means and variances through the function 
    \code{\link{priorControl}}.
  }
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
  \item{warn}{
    A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when
    estimation fails. 
    The default is controlled by \code{\link{mclust.options}}.
  }
  \item{verbose}{
    A logical controlling if a text progress bar is displayed during the
    fitting procedure. By default is \code{TRUE} if the session is 
    interactive, and \code{FALSE} otherwise.
  }  
  \item{\dots }{Further arguments passed to or from other methods.}
}

\value{
 An object of class \code{'MclustSSC'} providing the optimal (according 
 to BIC) Gaussian mixture model for semi-supervised classification.
 
 The details of the output components are as follows:

  \item{call}{The matched call.} 

  \item{data}{The input data matrix.} 
  
  \item{class}{The input class labels (including \code{NA}s for unknown labels.} 
    
  \item{modelName}{A character string specifying the "best" estimated model.}
    
  \item{G}{A numerical value specifying the number of mixture components or classes of the "best" estimated model.} 
    
  \item{n}{The total number of observations in the data.}

  \item{d}{The dimension of the data.}

  \item{BIC}{All BIC values.}

  \item{loglik}{Log-likelihood for the selected model.}
  
  \item{df}{Number of estimated parameters.}

  \item{bic}{Optimal BIC value.}
  
  \item{parameters}{
   A list with the following components:
   \describe{
    \item{\code{pro}}{
        A vector whose \emph{k}th component is the mixing proportion for
        the \emph{k}th component of the mixture model.
    }
    \item{\code{mean}}{
        The mean for each component. If there is more than one component,
        this is a matrix whose kth column is the mean of the \emph{k}th
        component of the mixture model.
    }
    \item{\code{variance}}{
        A list of variance parameters for the model.
        The components of this list depend on the model specification. 
        See the help file for \code{\link{mclustVariance}} for details.
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
  
  \item{prior}{
  The prior used (if any).
  }
  \item{control}{
  A list of control parameters used in the EM algorithm. 
  }
}

\details{
The semi-supervised approach implemented in \code{MclustSSC()} is a simple Gaussian mixture model for classification where at the first M-step only observations with known class labels are used for parameters estimation. Then, a standard EM algorithm is used for updating the probabiltiy of class membership for unlabelled data while keeping fixed the probabilities for labelled data. 
}

\references{
Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317. 
}

\author{Luca Scrucca}

\seealso{
  \code{\link{summary.MclustSSC}}, 
  \code{\link{plot.MclustSSC}}, 
  \code{\link{predict.MclustSSC}}
}

\examples{
# Simulate two overlapping groups
n <- 200
pars <- list(pro = c(0.5, 0.5),
             mean = matrix(c(-1,1), nrow = 2, ncol = 2, byrow = TRUE),
             variance = mclustVariance("EII", d = 2, G = 2))
pars$variance$sigmasq <- 1
data <- sim("EII", parameters = pars, n = n, seed = 12)
class <- data[,1]
X <- data[,-1]
clPairs(X, class, symbols = c(1,2), main = "Full classified data")

# Randomly remove labels
cl <- class; cl[sample(1:n, size = 195)] <- NA
table(cl, useNA = "ifany")
clPairs(X, ifelse(is.na(cl), 0, class),
        symbols = c(0, 16, 17), colors = c("grey", 4, 2),
        main = "Partially classified data")

# Fit semi-supervised classification model
mod_SSC  <- MclustSSC(X, cl)
summary(mod_SSC, parameters = TRUE)

pred_SSC <- predict(mod_SSC)
table(Predicted = pred_SSC$classification, Actual = class)

ngrid <- 50
xgrid <- seq(-3, 3, length.out = ngrid)
ygrid <- seq(-4, 4.5, length.out = ngrid)
xygrid <- expand.grid(xgrid, ygrid)
pred_SSC  <- predict(mod_SSC, newdata = xygrid)
col <- mclust.options("classPlotColors")[class]
pch <- class
pch[!is.na(cl)] = ifelse(cl[!is.na(cl)] == 1, 19, 17)
plot(X, pch = pch, col = col)
contour(xgrid, ygrid, matrix(pred_SSC$z[,1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lty = 2, lwd = 2)
}
\keyword{classification}
