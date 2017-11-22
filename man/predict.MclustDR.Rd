\name{predict.MclustDR}
\alias{predict.MclustDR}
\alias{predict2D.MclustDR}


\title{Classify multivariate observations on a dimension reduced subspace by Gaussian finite mixture modeling}

\description{Classify multivariate observations on a dimension reduced subspace estimated from a Gaussian finite mixture model.}

\usage{
  \method{predict}{MclustDR}(object, dim = 1:object$numdir, newdata, eval.points, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'MclustDR'} resulting from a call to \code{\link{MclustDR}.}}

  \item{dim}{the dimensions of the reduced subspace used for prediction.}

  \item{newdata}{a data frame or matrix giving the data. If missing the data obtained from the call to \code{\link{MclustDR}} are used.}

  \item{eval.points}{a data frame or matrix giving the data projected on the reduced subspace. If provided \code{newdata} is not used.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
Returns a list of with the following components:

  \item{dir}{a matrix containing the data projected onto the \code{dim} dimensions of the reduced subspace.}
  \item{density}{densities from mixture model for each data point.}
  \item{z}{a matrix whose \emph{[i,k]}th entry is the probability that 
           observation \emph{i} in \code{newdata} belongs to the \emph{k}th class.}
  \item{uncertainty}{The uncertainty associated with the classification.}
  \item{classification}{A vector of values giving the MAP classification.}
}

\references{
  Scrucca, L. (2010) Dimension reduction for model-based clustering. 
  \emph{Statistics and Computing}, 20(4), pp. 471-484.
}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{MclustDR}}.}

\examples{
mod = Mclust(iris[,1:4])
dr = MclustDR(mod)
pred = predict(dr)
str(pred)

data(banknote)
mod = MclustDA(banknote[,2:7], banknote$Status)
dr = MclustDR(mod)
pred = predict(dr)
str(pred)
}

\keyword{multivariate}
