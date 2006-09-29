\name{cross}
\alias{cross}
\non_function{}
\title{Simulated Cross Data}
\usage{data(cross)}
\description{
A 500 by 3 matrix in which the first column is the classification and
the remaining columns are two data from a simulation of two crossed
elliptical Gaussians. 
}

\examples{
# This dataset was created as follows
\dontrun{
n <- 250 
set.seed(0)
cross <- rbind(matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9)),
           matrix(rnorm(n*2), n, 2) \%*\% diag(c(1,9))[,2:1])
cross <- cbind(c(rep(1,n),rep(2,n)), x)
}
}

\keyword{datasets}
