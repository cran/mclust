\name{plot.MclustSSC}
\alias{plot.MclustSSC}

\title{Plotting method for MclustSSC semi-supervised classification}

\description{
Plots for semi-supervised classification based on Gaussian finite mixture models.
}

\usage{
\method{plot}{MclustSSC}(x, what = c("BIC", "classification", "uncertainty"), \dots)
}

\arguments{
  \item{x}{
    An object of class \code{'MclustSSC'} resulting from a call to \code{\link{MclustSSC}}.
  }

 \item{what}{
    A string specifying the type of graph requested. Available choices are:
    \describe{
      \item{\code{"BIC"} =}{plot of BIC values used for model selection, i.e. for choosing the model class covariances.}
      \item{\code{"classification"} =}{a plot of data with points marked based on the known and the predicted classification.}
      \item{\code{"uncertainty"} =}{a plot of classification uncertainty.}
   }
   If not specified, in interactive sessions a menu of choices is proposed.
  }

  \item{\dots}{further arguments passed to or from other methods. See  \code{\link{plot.Mclust}}.}
}

%\value{}
  
%\details{}

\author{Luca Scrucca}

\seealso{
  \code{\link{MclustSSC}}
}

\examples{
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
summary(mod_SSC, parameters = TRUE)

pred_SSC <- predict(mod_SSC)
table(Predicted = pred_SSC$classification, Actual = class, useNA = "ifany")

plot(mod_SSC, what = "BIC")
plot(mod_SSC, what = "classification")
plot(mod_SSC, what = "uncertainty")
}

\keyword{multivariate}
