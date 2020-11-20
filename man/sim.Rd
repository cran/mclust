\name{sim}
\alias{sim}
\title{
  Simulate from Parameterized MVN Mixture Models
}
\description{
  Simulate data from parameterized MVN mixture models.
}
\usage{
sim(modelName, parameters, n, seed = NULL, \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model. The help file for
    \code{\link{mclustModelNames}} describes the available models.
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
  \item{n}{
      An integer specifying the number of data points to be simulated.
  }
  \item{seed}{
    An optional integer argument to \code{set.seed} for reproducible
    random class assignment.  By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{sim}.
  }
  \item{\dots }{
     Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}
\value{
  A matrix in which first column is the classification and the remaining
  columns are the \code{n} observations simulated from the specified MVN 
  mixture model.
\item{Attributes:}{
      \code{"modelName"} A character string indicating the variance 
      model used for the simulation.
  }
}
\details{
  This function can be used with an indirect or list call using
  \code{do.call}, allowing the output of e.g. \code{mstep}, \code{em},
  \code{me}, \code{Mclust} to be passed directly without the need to
  specify individual parameters as arguments. 
}

\seealso{
  \code{\link{simE}}, \dots,
  \code{\link{simVVV}},
  \code{\link{Mclust}},
  \code{\link{mstep}},
  \code{\link{do.call}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
irisModel <- mclustModel(iris[,-5], irisBIC)
names(irisModel)
irisSim <- sim(modelName = irisModel$modelName, 
               parameters = irisModel$parameters, 
               n = nrow(iris))

\dontrun{
  do.call("sim", irisModel) # alternative call
}

par(pty = "s", mfrow = c(1,2))

dimnames(irisSim) <- list(NULL, c("dummy", (dimnames(iris)[[2]])[-5]))

dimens <- c(1,2)
lim1 <- apply(iris[,dimens],2,range)
lim2 <- apply(irisSim[,dimens+1],2,range)
lims <- apply(rbind(lim1,lim2),2,range)
xlim <- lims[,1]
ylim <- lims[,2]

coordProj(iris[,-5], parameters=irisModel$parameters, 
          classification=map(irisModel$z), 
          dimens=dimens, xlim=xlim, ylim=ylim)

coordProj(iris[,-5], parameters=irisModel$parameters, 
          classification=map(irisModel$z), truth = irisSim[,-1],
          dimens=dimens, xlim=xlim, ylim=ylim)

irisModel3 <- mclustModel(iris[,-5], irisBIC, G=3)
irisSim3 <- sim(modelName = irisModel3$modelName, 
               parameters = irisModel3$parameters, n = 500, seed = 1)
\dontrun{
 irisModel3$n <- NULL
 irisSim3 <- do.call("sim",c(list(n=500,seed=1),irisModel3)) # alternative call
}
clPairs(irisSim3[,-1], cl = irisSim3[,1])
}
\keyword{cluster}
