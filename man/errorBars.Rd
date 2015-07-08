\name{errorBars}
\alias{errorBars}

\title{Draw error bars on a plot}

\description{
Draw error bars at x from upper to lower. If \code{horizontal = FALSE} (default)
bars are drawn vertically, otherwise horizontally.
}

\usage{
errorBars(x, upper, lower, width = 0.1, code = 3, angle = 90, horizontal = FALSE, \dots)
}

\arguments{
  \item{x}{A vector of values where the bars must be drawn.}
  \item{upper}{A vector of upper values where the bars must end.}
  \item{lower}{A vector of lower values where the bars must start.}
  \item{width}{A value specifying the width of the end-point segment.}
  \item{code}{An integer code specifying the kind of arrows to be drawn. For details see \code{\link[graphics]{arrows}}.}
  \item{angle}{A value specifying the angle at the arrow edge. For details see \code{\link[graphics]{arrows}}.}
  \item{horizontal}{A logical specifying if bars should be drawn vertically (default) or horizontally.}
  \item{\dots}{Further arguments are passed to \code{\link[graphics]{arrows}}.}
}

%\value{}

\examples{
par(mfrow=c(2,2))
# Create a simple example dataset
x <- 1:5
n <- c(10, 15, 12, 6, 3)
se <- c(1, 1.2, 2, 1, .5)
# upper and lower bars
b <- barplot(n, ylim = c(0, max(n)*1.5))
errorBars(b, lower = n-se, upper = n+se, lwd = 2, col = "red3")
# one side bars
b <- barplot(n, ylim = c(0, max(n)*1.5))
errorBars(b, lower = n, upper = n+se, lwd = 2, col = "red3", code = 1)
# 
plot(x, n, ylim = c(0, max(n)*1.5), pch = 0)
errorBars(x, lower = n-se, upper = n+se, lwd = 2, col = "red3")
#
dotchart(n, labels = x, pch = 19, xlim = c(0, max(n)*1.5))
errorBars(x, lower = n-se, upper = n+se, col = "red3", horizontal = TRUE)
}
