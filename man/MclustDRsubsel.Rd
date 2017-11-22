\name{MclustDRsubsel}
\alias{MclustDRsubsel}
\alias{print.MclustDRsubsel}
\alias{MclustDRsubsel_classif}
\alias{MclustDRsubsel_cluster}
\alias{MclustDRrecoverdir}
\alias{MclustDRsubsel1cycle}
\alias{print.MclustDRsubsel}
\alias{summary.MclustDRsubsel}

\title{Subset selection for GMMDR directions based on BIC.}

\description{
Implements a subset selection method for selecting the relevant directions spanning the dimension reduction subspace for visualizing the clustering or classification structure obtained from a finite mixture of Gaussian densities.}

\usage{
MclustDRsubsel(object, G = 1:9,
                       modelNames = mclust.options("emModelNames"), 
                       \dots,
                       bic.stop = 0, bic.cutoff = 0, 
                       mindir = 1, 
                       verbose = interactive())
}

\arguments{
  \item{object}{An object of class \code{'MclustDR'} resulting from a call to \code{\link{MclustDR}}.}
  \item{G}{An integer vector specifying the numbers of mixture components or clusters.}
  \item{modelNames}{A vector of character strings indicating the models to be fitted. See \code{\link{mclustModelNames}} for a description of the available models.}
  \item{\dots}{Further arguments passed through \code{\link{Mclust}} or \code{\link{MclustDA}}.}
  \item{bic.stop}{A criterion to terminate the search. If maximal BIC difference is less than \code{bic.stop} then the algorithm stops. \cr
  Two tipical values are:
  \tabular{lcl}{
  \code{0}    \tab = \tab algorithm stops when the BIC difference becomes negative (default) \cr
  \code{-Inf} \tab = \tab algorithm continues until all directions have been selected }
  }
  \item{bic.cutoff}{A value specifying how to select simplest ``best'' model within \code{bic.cutoff} from the maximum value achieved. Setting this to \code{0} (default) simply select the model with the largest BIC difference.}
  \item{mindir}{An integer value specifying the minimum number of directions to be estimated.}
  \item{verbose}{A logical or integer value specifying if and how much detailed information should be reported during the iterations of the algorithm. \cr
  Possible values are:
  \tabular{ll}{  
  \code{0} or \code{FALSE} \tab = no trace info is shown; \cr
  \code{1} or \code{TRUE} \tab = a trace info is shown at each step of the search; \cr
  \code{2} \tab = a more detailed trace info is is shown.}
  }
}
\details{
The GMMDR method aims at reducing the dimensionality by identifying a set of linear combinations, ordered by importance as quantified by the associated eigenvalues, of the original features which capture most of the clustering or classification structure contained in the data. This is implemented in \code{\link{MclustDR}}.

The \code{MclustDRsubsel} function implements the greedy forward search algorithm discussed in Scrucca (2010) to prune the set of all GMMDR directions. The criterion used to select the relevant directions is based on the BIC difference between a clustering model and a model in which the feature proposal has no clustering relevance. The steps are the following:

1. Select the first feature to be the one which maximizes the BIC difference between the best clustering model and the model which assumes no clustering, i.e. a single component.

2. Select the next feature amongst those not previously included, to be the one which maximizes the BIC difference.

3. Iterate the previous step until all the BIC differences for the inclusion of a feature become less than \code{bic.stop}.

At each step, the search over the model space is performed with respect to the model parametrisation and the number of clusters. 
}

\value{
  An object of class \code{'MclustDRsubsel'} which inherits from \code{'MclustDR'}, so it has the same components of the latter plus the following:

  \item{basisx}{The basis of the estimated dimension reduction subspace expressed in terms of the original variables.}

  \item{std.basisx}{The basis of the estimated dimension reduction subspace expressed in terms of the original variables standardized to have unit standard deviation.}

}

\references{
Scrucca, L. (2010) Dimension reduction for model-based clustering. \emph{Statistics and Computing}, 20(4), pp. 471-484.

Scrucca, L. (2014) Graphical Tools for Model-based Mixture Discriminant Analysis. \emph{Advances in Data Analysis and Classification}, 8(2), pp. 147-165
}

\author{Luca Scrucca}
%\note{}

\seealso{
\code{\link{MclustDR}}, \code{\link{Mclust}}, \code{\link{MclustDA}}.
}

\examples{
\dontrun{
# clustering
data(crabs, package = "MASS")
x = crabs[,4:8]
class = paste(crabs$sp, crabs$sex, sep = "|")
mod = Mclust(x)
table(class, mod$classification)
dr = MclustDR(mod)
summary(dr)
plot(dr)
drs = MclustDRsubsel(dr)
summary(drs)
table(class, drs$class)
plot(drs, what = "scatterplot")
plot(drs, what = "pairs")
plot(drs, what = "contour")
plot(drs, what = "boundaries")
plot(drs, what = "evalues")

# classification
data(banknote)
da = MclustDA(banknote[,2:7], banknote$Status)
table(banknote$Status, predict(da)$class)
dr = MclustDR(da)
summary(dr)
drs = MclustDRsubsel(dr)
summary(drs)
table(banknote$Status, predict(drs)$class)
plot(drs, what = "scatterplot")
plot(drs, what = "classification")
plot(drs, what = "boundaries")}
}

\keyword{multivariate}