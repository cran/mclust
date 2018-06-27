\name{hdrlevels}
\alias{hdrlevels}

\title{Highest Density Region (HDR) Levels}

\description{
Compute the levels of Highest Density Regions (HDRs) for any density and probability levels.
}

\usage{
hdrlevels(density, prob)
}

\arguments{
  \item{density}{A vector of density values computed on a set of (observed) evaluation points.}
  \item{prob}{A vector of probability levels in the range \eqn{[0,1]}.}
}

\value{
The function returns a vector of density values corresponding to HDRs at given probability levels.
} 

\details{
From Hyndman (1996), let \eqn{f(x)} be the density function of a random 
variable \eqn{X}. Then the \eqn{100(1-\alpha)\%} HDR is the subset 
\eqn{R(f_\alpha)} of the sample space of \eqn{X} such that
\deqn{
R(f_\alpha) = {x : f(x) \ge f_\alpha }
}
where \eqn{f_\alpha} is the largest constant such that 
\eqn{
Pr( X \in R(f_\alpha)) \ge 1-\alpha
}
}

\seealso{
\code{\link{plot.densityMclust}}
}

\references{
Rob J. Hyndman (1996) Computing and Graphing Highest Density Regions. \emph{The American Statistician}, 50(2):120-126.
}

\author{L. Scrucca}

\examples{
# Example: univariate Gaussian
x <- rnorm(1000)
f <- dnorm(x)
a <- c(0.5, 0.25, 0.1)
(f_a <- hdrlevels(f, prob = 1-a))

plot(x, f)
abline(h = f_a, lty = 2)
text(max(x), f_a, labels = paste0("f_", a), pos = 3)

mean(f > f_a[1])
range(x[which(f > f_a[1])])
qnorm(1-a[1]/2)

mean(f > f_a[2])
range(x[which(f > f_a[2])])
qnorm(1-a[2]/2)

mean(f > f_a[3])
range(x[which(f > f_a[3])])
qnorm(1-a[3]/2)

# Example 2: univariate Gaussian mixture
set.seed(1)
cl <- sample(1:2, size = 1000, prob = c(0.7, 0.3), replace = TRUE)
x <- ifelse(cl == 1, 
            rnorm(1000, mean = 0, sd = 1),
            rnorm(1000, mean = 4, sd = 1))
f <- 0.7*dnorm(x, mean = 0, sd = 1) + 0.3*dnorm(x, mean = 4, sd = 1)

a <- 0.25
(f_a <- hdrlevels(f, prob = 1-a))

plot(x, f)
abline(h = f_a, lty = 2)
text(max(x), f_a, labels = paste0("f_", a), pos = 3)

mean(f > f_a)

# find the regions of HDR
ord <- order(x)
f <- f[ord]
x <- x[ord]
x_a <- x[f > f_a]
j <- which.max(diff(x_a))
region1 <- x_a[c(1,j)]
region2 <- x_a[c(j+1,length(x_a))]
plot(x, f, type = "l")
abline(h = f_a, lty = 2)
abline(v = region1, lty = 3, col = 2)
abline(v = region2, lty = 3, col = 3)
}
\keyword{density}
