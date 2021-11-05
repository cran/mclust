\name{predict.MclustDA}
\alias{predict.MclustDA}

\title{Classify multivariate observations by Gaussian finite mixture modeling}

\description{Classify multivariate observations based on Gaussian finite mixture models estimated by \code{\link{MclustDA}}.}

\usage{
  \method{predict}{MclustDA}(object, newdata, prop = object$prop, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.}

  \item{newdata}{a data frame or matrix giving the data. If missing the train data obtained from the call to \code{\link{MclustDA}} are classified.}

  \item{prop}{the class proportions or prior class probabilities to belong to each class; by default, this is set at the class proportions in the training data.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
Returns a list of with the following components:
  \item{classification}{a factor of predicted class labels for \code{newdata}.}
  \item{z}{a matrix whose \emph{[i,k]}th entry is the probability that 
           observation \emph{i} in \code{newdata} belongs to the \emph{k}th class.}
}

\author{Luca Scrucca}

% \note{}

\seealso{\code{\link{MclustDA}}.}

\examples{
\donttest{
odd <- seq(from = 1, to = nrow(iris), by = 2)
even <- odd + 1
X.train <- iris[odd,-5]
Class.train <- iris[odd,5]
X.test <- iris[even,-5]
Class.test <- iris[even,5]

irisMclustDA <- MclustDA(X.train, Class.train)

predTrain <- predict(irisMclustDA)
predTrain
predTest <- predict(irisMclustDA, X.test)
predTest
}
}

\keyword{multivariate}
