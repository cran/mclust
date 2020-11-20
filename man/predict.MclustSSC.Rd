\name{predict.MclustSSC}
\alias{predict.MclustSSC}

\title{Classification of multivariate observations by semi-supervised Gaussian finite mixtures}

\description{Classify multivariate observations based on Gaussian finite mixture models estimated by \code{\link{MclustSSC}}.}

\usage{
  \method{predict}{MclustSSC}(object, newdata, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'MclustSSC'} resulting from a call to \code{\link{MclustSSC}}.}

  \item{newdata}{a data frame or matrix giving the data. If missing the train data obtained from the call to \code{\link{MclustSSC}} are classified.}

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

\seealso{\code{\link{MclustSSC}}.}

\examples{
\dontrun{
X <- iris[,1:4]
class <- iris$Species
# randomly remove class labels
set.seed(123)
class[sample(1:length(class), size = 120)] <- NA
table(class, useNA = "ifany")
clPairs(X, ifelse(is.na(class), 0, class),
        symbols = c(0, 16, 17, 18), colors = c("grey", 4, 2, 3),
        main = "Partially classified data")

# Fit semi-supervised classification model
mod_SSC  <- MclustSSC(X, class)

pred_SSC <- predict(mod_SSC)
table(Predicted = pred_SSC$classification, Actual = class, useNA = "ifany")

X_new = data.frame(Sepal.Length = c(5, 8),
                   Sepal.Width  = c(3.1, 4),
                   Petal.Length = c(2, 5),
                   Petal.Width  = c(0.5, 2))
predict(mod_SSC, newdata = X_new)
}
}

\keyword{classification}
