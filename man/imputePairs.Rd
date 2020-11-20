\name{imputePairs}
\alias{imputePairs}
\title{
  Pairwise Scatter Plots showing Missing Data Imputations
}
\description{
  Creates a scatter plot for each pair of variables in given data,
  allowing display of imputations for missing values in different
  colors and symbols than non missing values.
}
\usage{
imputePairs(data, dataImp, 
            symbols = c(1,16), colors = c("black", "red"), labels,
            panel = points, \dots, lower.panel = panel, upper.panel = panel, 
            diag.panel = NULL, text.panel = textPanel, label.pos = 0.5 + 
            has.diag/3, cex.labels = NULL, font.labels = 1, row1attop = TRUE, 
            gap = 0.2) 
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations containing
    missing values. Categorical variables are not allowed. If a matrix
    or data frame, rows correspond to observations and columns
    correspond to variables. 
  }
  \item{dataImp}{
    The dataset \code{data} with missing values imputed.
  }
  \item{symbols}{
    Either an integer or character vector assigning plotting symbols to
    the nonmissing data and impued  values, respectively. The default is a
    closed circle for the nonmissing data
    and an open circle for the imputed values.
  }
  \item{colors}{
    Either an integer or character vector assigning colors to
    the nonmissing data and impued  values, respectively. The default is 
    black for the nonmissing data and red for the imputed values.
  }
  \item{labels}{
    As in function \code{pairs}.
  }
  \item{panel}{
    As in function \code{pairs}.
  }
  \item{\dots}{
    As in function \code{pairs}.
  }
  \item{lower.panel}{
    As in function \code{pairs}.
  }
  \item{upper.panel}{
    As in function \code{pairs}.
  }
  \item{diag.panel}{
    As in function \code{pairs}.
  }
  \item{text.panel}{
    As in function \code{pairs}.
  }
  \item{label.pos}{
    As in function \code{pairs}.
  }
  \item{cex.labels}{
    As in function \code{pairs}.
  }
  \item{font.labels}{
    As in function \code{pairs}.
  }
  \item{row1attop}{
    As in function \code{pairs}.
  }
  \item{gap}{
    As in function \code{pairs}.
  }
}
\value{
  A pairs plot displaying the location of missing and nonmissing values.
}

\references{
  Schafer J. L. (1997). Analysis of Imcomplete Multivariate Data, Chapman and Hall.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{imputeData}}
}
\examples{
\dontrun{
# Note that package 'mix' must be installed
data(stlouis, package = "mix")
 
# impute the continuos variables in the stlouis data
stlimp <- imputeData(stlouis[,-(1:3)])

# plot imputed values
imputePairs(stlouis[,-(1:3)], stlimp)
}
}
\keyword{cluster}
