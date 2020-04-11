\name{dmvnorm}
\alias{dmvnorm}

\title{Density of multivariate Gaussian distribution}

\description{
Efficiently computes the density of observations for a generic multivariate Gaussian distribution.
}

\usage{
dmvnorm(data, mean, sigma, log = FALSE)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
 \item{mean}{
    A vector of means for each variable.
  }
 \item{sigma}{
    A positive definite covariance matrix.
  }
  \item{log}{
    A logical value indicating whether or not the logarithm of the densities 
    should be returned.
  }
}
\value{
  A numeric vector whose \emph{i}th element gives the density of the
  \emph{ith} observation in \code{data} for the multivariate Gaussian 
	distribution with parameters \code{mean} and \code{sigma}.
}

\seealso{
  \code{\link{dnorm}},
  \code{\link{dens}}
}

\examples{
# univariate
ngrid <- 101
x <- seq(-5, 5, length = ngrid)
dens <- dmvnorm(x, mean = 1, sigma = 5)
plot(x, dens, type = "l")

# bivariate
ngrid <- 101
x1 <- x2 <- seq(-5, 5, length = ngrid)
mu <- c(1,0)
sigma <- matrix(c(1,0.5,0.5,2), 2, 2)
dens <- dmvnorm(as.matrix(expand.grid(x1, x2)), mu, sigma)
dens <- matrix(dens, ngrid, ngrid)
image(x1, x2, dens)
contour(x1, x2, dens, add = TRUE)
}

