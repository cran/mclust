\name{entPlot}
\alias{entPlot}
\title{
Plot Entropy Plots
}
\description{
Plot "entropy plots" to help select the number of classes from a hierarchy of combined clusterings.
}
\usage{
entPlot(z, combiM, abc = c("standard", "normalized"), reg = 2, ...)
}
\arguments{
  \item{z}{
	A matrix whose \code{[i,k]}th entry is the probability that observation \emph{i} in the data belongs to the \emph{k}th class, for the initial solution (ie before any combining). Typically, the one returned by \code{Mclust}/BIC.
	}
  \item{combiM}{
	A list of "combining matrices" (as provided by \code{clustCombi}), ie \code{combiM[[K]]} is the matrix whose \emph{k}th row contains only zeros, but in columns corresponding to the labels of the classes in the \emph{(K+1)}-classes solution to be merged to get the \emph{K}-classes combined solution. \code{combiM} must contain matrices from \code{K} = number of classes in \code{z} to one.
	}
  \item{abc}{
	Choose one or more of: "standard", "normalized", to specify whether the number of observations involved in each combining step should be taken into account to scale the plots or not.
	}
  \item{reg}{
	The number of parts of the piecewise linear regression for the entropy plots. Choose one or more of: 2 (for 1 change-point), 3 (for 2 change-points).
	}
  \item{\dots}{
	Other graphical arguments to be passed to the plot functions.
	}
}
\details{
Please see the article cited in the references for more details. A clear elbow in the "entropy plot" should suggest the user to consider the corresponding number(s) of class(es).
}
\value{
 if \code{abc = "standard"}, plots the entropy against the number of clusters and the difference between the entropy of successive combined solutions against the number of clusters.
 if \code{abc = "normalized"}, plots the entropy against the cumulated number of observations involved in the successive combining steps and the difference between the entropy of successive combined solutions divided by the number of observations involved in the corresponding combining step against the number of clusters.
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{plot.clustCombi}}, \code{\link{combiPlot}}, \code{\link{clustCombi}}
}
\examples{
\dontrun{
data(Baudry_etal_2010_JCGS_examples)
# run Mclust to get the MclustOutput
output <- clustCombi(ex4.2, modelNames = "VII") 

entPlot(output$MclustOutput$z, output$combiM, reg = c(2,3)) 
# legend: in red, the single-change-point piecewise linear regression;
#         in blue, the two-change-point piecewise linear regression.
}
}
\keyword{ cluster }
