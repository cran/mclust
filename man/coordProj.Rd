\name{coordProj}
\alias{coordProj}
\title{
  Coordinate projections of data in more than two dimensions modelled by an MVN
  mixture.
}
\description{
  Plots coordinate projections given data in more than two dimensions
  and parameters of an MVN mixture model for the data.
}
\usage{
coordProj(data, \dots, dimens = c(1, 2),
          type = c("classification","uncertainty","errors"), ask = TRUE,
          quantiles = c(0.75, 0.95), symbols, scale = FALSE,
          identify = FALSE, CEX = 1, PCH = ".", xlim, ylim)
}
\arguments{
  \item{data}{
    A numeric matrix or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables.
  }
  \item{dimens}{
    A vector of length 2 giving the integer dimensions of the
    desired coordinate projections. The default is
    \code{c(1,2)}, in which the first
    dimension is plotted against the second.
  }
  \item{\dots}{
    One or more of the following:
    \describe{
      \item{classification}{
	A numeric or character vector representing a classification of
	observations (rows) of \code{data}.
      }
      \item{uncertainty}{
	A numeric vector of values in \emph{(0,1)} giving the
	uncertainty of each data point.
      }
      \item{z}{
	A matrix in which the \code{[i,k]}th entry gives the
	probability of observation \emph{i} belonging to the \emph{k}th class. 
	Used to compute \code{classification} and
	\code{uncertainty} if those arguments aren't available.
      }
      \item{truth}{
	A numeric or character vector giving a known
	classification of each data point.
	If \code{classification}
	or \code{z} is also present, 
	this is used for displaying classification errors.
      }
      \item{mu}{
	A matrix whose columns are the means of each group. 
      }
      \item{sigma}{
	A three dimensional array 
	in which \code{sigma[,,k]} gives the covariance
	for the \emph{k}th group.
      }
      \item{decomp}{
	A list with \code{scale}, \code{shape} and \code{orientation}
	components giving an alternative form for the covariance
	structure of the mixture model. 
      }
    }
  }
  \item{type}{
    Any subset of 
    \code{c("classification","uncertainty","errors")}. The function will
    produce the corresponding plot if it has been supplied sufficient
    information to do so. If more than one plot is possible then users
    will be asked to choose from a menu if \code{ask=TRUE}. 
  }
  \item{ask}{
    A logical variable indicating whether or not a menu should be
    produced when more than one plot is possible. The default is
    \code{ask=TRUE}. 
  }
  \item{quantiles}{
    A vector of length 2 giving quantiles used in plotting
    uncertainty. The smallest symbols correspond to the smallest
    quantile (lowest uncertainty), medium-sized (open) symbols to points
    falling between the given quantiles, and large (filled) symbols to
    those in the largest quantile (highest uncertainty). The default is
    \emph{(0.75,0.95)}. 
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to
    each unique class in \code{classification}. Elements in
    \code{symbols} correspond to classes in \code{classification} in
    sorted order. Default: If \emph{G} is the number of groups in the
    classification, the first \emph{G} symbols in
    \code{.Mclust\$symbols}, otherwise if \emph{G} is less than 27 then
    the first \emph{G} capital letters in the Roman alphabet. 
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=FALSE} 
  }
  \item{identify}{
    A logical variable indicating whether or not to add a title to the plot
    identifying the dimensions used.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{PCH}{
    An argument specifying the symbol to be used when a classificatiion
    has not been specified for the data. The default value is a small dot ".".
  }
  \item{xlim, ylim}{
    Arguments specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
}
\section{Side Effects}{
  Coordinate projections of the data, possibly showing location of the
  mixture components, classification, uncertainty, and/or classification
  errors. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}.
  
  
  C. Fraley and A. E. Raftery (2002).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis.
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{clPairs}},
  \code{\link{randProj}},
  \code{\link{mclust2Dplot}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

msEst <- mstepVVV(irisMatrix, unmap(irisClass))

par(pty = "s", mfrow = c(1,2))
coordProj(irisMatrix,dimens=c(2,3), truth = irisClass, 
          mu = msEst$mu, sigma = msEst$sigma, z = msEst$z)
do.call("coordProj", c(list(data=irisMatrix, dimens=c(2,3), truth=irisClass),
                       msEst))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
