\name{MclustDR}
\alias{MclustDR}
\alias{print.MclustDR}

\title{Dimension reduction for model-based clustering and classification}

\description{
A dimension reduction method for visualizing the clustering or classification structure obtained from a finite mixture of Gaussian densities.
}

\usage{
MclustDR(object, lambda = 1, normalized = TRUE, Sigma,
         tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{object}{An object of class \code{'Mclust'} or \code{'MclustDA'} 
  resulting from a call to, respectively, \code{\link{Mclust}} or 
  \code{\link{MclustDA}}.}
  
  \item{lambda}{A tuning parameter in the range [0,1] as described in 
  Scrucca (2014). The directions that mostly separate the estimated clusters 
  or classes are recovered using the default value 1. Users can set this 
  parameter to balance the relative importance of information derived from 
  cluster/class means and covariances. For instance, a value of 0.5 gives 
  equal importance to differences in means and covariances among clusters/classes.}
  
  \item{normalized}{Logical. If \code{TRUE} directions are normalized to unit norm.}
  
  \item{Sigma}{Marginal covariance matrix of data. If not provided is estimated by the MLE of observed data.}
  
  \item{tol}{A tolerance value.}
}

\details{
The method aims at reducing the dimensionality by identifying a set of linear combinations, ordered by importance as quantified by the associated eigenvalues, of the original features which capture most of the clustering or classification structure contained in the data. 

Information on the dimension reduction subspace is obtained from the variation on group means and, depending on the estimated mixture model, on the variation on group covariances (see Scrucca, 2010). 

Observations may then be projected onto such a reduced subspace, thus providing summary plots which help to visualize the underlying structure.

The method has been extended to the supervised case, i.e. when the true classification is known (see Scrucca, 2014).

This implementation doesn't provide a formal procedure for the selection of dimensionality. A future release will include one or more methods.
}

\value{
  An object of class \code{'MclustDR'} with the following components:

  \item{call}{The matched call} 

  \item{type}{A character string specifying the type of model for which the dimension reduction is computed. Currently, possible values are \code{"Mclust"} for clustering, and \code{"MclustDA"} or \code{"EDDA"} for classification.} 
    
  \item{x}{The data matrix.}

  \item{Sigma}{The covariance matrix of the data.}

  \item{mixcomp}{A numeric vector specifying the mixture component of each data observation.}
  
  \item{class}{A factor specifying the classification of each data observation. For model-based clustering this is equivalent to the corresponding mixture component. For model-based classification this is the known classification.}
    
  \item{G}{The number of mixture components.}
  
  \item{modelName}{The name of the parameterization of the estimated mixture model(s). See \code{\link{mclustModelNames}}.}

  \item{mu}{A matrix of means for each mixture component.}
    
  \item{sigma}{An array of covariance matrices for each mixture component.}
  
  \item{pro}{The estimated prior for each mixture component.}

  \item{M}{The kernel matrix.}
  
  \item{lambda}{The tuning parameter.}
  
  \item{evalues}{The eigenvalues from the generalized eigen-decomposition of the kernel matrix.}

  \item{raw.evectors}{The raw eigenvectors from the generalized eigen-decomposition of the kernel matrix, ordered according to the eigenvalues.}

  \item{basis}{The basis of the estimated dimension reduction subspace.}

  \item{std.basis}{The basis of the estimated dimension reduction subspace standardized to variables having unit standard deviation.}

  \item{numdir}{The dimension of the projection subspace.}

  \item{dir}{The estimated directions, i.e. the data projected onto the estimated dimension reduction subspace.}

}

\references{
Scrucca, L. (2010) Dimension reduction for model-based clustering. \emph{Statistics and Computing}, 20(4), pp. 471-484.

Scrucca, L. (2014) Graphical Tools for Model-based Mixture Discriminant Analysis. \emph{Advances in Data Analysis and Classification}, 8(2), pp. 147-165.
}

\author{Luca Scrucca}
%\note{}

\seealso{
\code{\link{summary.MclustDR}}, \code{\link{plot.MclustDR}}, \code{\link{Mclust}}, \code{\link{MclustDA}}.
}

\examples{
# clustering
data(diabetes)
mod <- Mclust(diabetes[,-1])
summary(mod)

dr <- MclustDR(mod)
summary(dr)
plot(dr, what = "scatterplot")
plot(dr, what = "evalues")

dr <- MclustDR(mod, lambda = 0.5) 
summary(dr)
plot(dr, what = "scatterplot")
plot(dr, what = "evalues")

# classification
data(banknote)

da <- MclustDA(banknote[,2:7], banknote$Status, modelType = "EDDA")
dr <- MclustDR(da)
summary(dr)

da <- MclustDA(banknote[,2:7], banknote$Status)
dr <- MclustDR(da)
summary(dr)
}

\keyword{multivariate}
