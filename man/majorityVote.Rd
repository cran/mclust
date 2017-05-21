\name{majorityVote}
\alias{majorityVote}

\title{Majority vote}

\description{
A function to compute the majority vote (some would say plurality) label in a vector of labels, breaking ties at random.}

\usage{
majorityVote(x)
}

\arguments{
  \item{x}{A vector of values, either numerical or not.}
}

\value{A list with the following components:
  \item{table}{A table of votes for each unique value of \code{x}.}
  \item{ind}{An integer specifying which unique value of \code{x} corresponds to the majority vote.}
  \item{majority}{A string specifying the majority vote label.}
}

%\seealso{}
\author{L. Scrucca}

\examples{
x <- c("A", "C", "A", "B", "C", "B", "A")
majorityVote(x)
}
