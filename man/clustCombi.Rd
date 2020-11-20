\name{clustCombi}
\alias{clustCombi}
\alias{print.clustCombi}
\alias{summary.clustCombi}
\alias{print.summary.clustCombi}

\title{
Combining Gaussian Mixture Components for Clustering
}
\description{
Provides a hierarchy of combined clusterings from the EM/BIC Gaussian mixture solution to one class, following the methodology proposed in the article cited in the references.
}
\usage{
clustCombi(object = NULL, data = NULL, \dots)
}
\arguments{
  \item{object}{
    An object returned by \code{\link{Mclust}} giving the optimal (according to BIC) parameters, conditional probabilities, and log-likelihood, together with the associated classification and its uncertainty. If not provided, the \code{data} argument must be specified.
  }
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables. If the \code{object} argument is not provided, the function \code{\link{Mclust}} is applied to the given \code{data} to fit a mixture model.}
  \item{\dots}{
    Optional arguments to be passed to called functions. Notably, any argument (such as the numbers of components for which the BIC is computed; the models to be fitted by EM; initialization parameters for the EM algorithm, etc.) to be passed to \code{\link{Mclust}} in case \code{object = NULL}. Please see the \code{\link{Mclust}} documentation for more details.
  }
}
\details{
Mclust provides a Gaussian mixture fitted to the data by maximum likelihood through the EM algorithm, for the model and number of components selected according to BIC. The corresponding components are hierarchically combined according to an entropy criterion, following the methodology described in the article cited in the references section. The solutions with numbers of classes between the one selected by BIC and one are returned as a \code{clustCombi} class object.
}
\value{
A list of class \code{clustCombi} giving the hierarchy of combined solutions from the number of components selected by BIC to one. The details of the output components are as follows:
 \item{classification}{A list of the data classifications obtained for each combined solution of the hierarchy through a MAP assignment}
 \item{combiM}{A list of matrices. \code{combiM[[K]]} is the matrix used to combine the components of the (K+1)-classes solution to get the K-classes solution. Please see the examples.}
 \item{combiz}{A list of matrices. \code{combiz[[K]]} is a matrix whose [i,k]th entry is the probability that observation i in the data belongs to the kth class according to the K-classes combined solution.}
 \item{MclustOutput}{A list of class \code{Mclust}. Output of a call to the Mclust function (as provided by the user or the result of a call to the Mclust function) used to initiate the combined solutions hierarchy: please see the \code{\link{Mclust}} function documentation for details.}
}
\references{
J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{plot.clustCombi}}
}
\examples{
data(Baudry_etal_2010_JCGS_examples)

# run Mclust using provided data
output <- clustCombi(data = ex4.1) 
\dontrun{
# or run Mclust and then clustcombi on the returned object
mod <- Mclust(ex4.1)
output <- clustCombi(mod)
}

output
summary(output)

\dontrun{
# run Mclust using provided data and any further optional argument provided
output <- clustCombi(data = ex4.1, modelName = "EEV", G = 1:15)
}

# plot the hierarchy of combined solutions
plot(output, what = "classification") 
# plot some "entropy plots" which may help one to select the number of classes
plot(output, what = "entropy") 
# plot the tree structure obtained from combining mixture components
plot(output, what = "tree") 

# the selected model and number of components obtained from Mclust using BIC
output$MclustOutput 

# the matrix whose [i,k]th entry is the probability that i-th observation in 
# the data belongs to the k-th class according to the BIC solution
head( output$combiz[[output$MclustOutput$G]] ) 
# the matrix whose [i,k]th entry is the probability that i-th observation in 
# the data belongs to the k-th class according to the first combined solution
head( output$combiz[[output$MclustOutput$G-1]] ) 
# the matrix describing how to merge the 6-classes solution to get the 
# 5-classes solution
output$combiM[[5]] 
# for example the following code returns the label of the class (in the 
# 5-classes combined solution) to which the 4th class (in the 6-classes
# solution) is assigned. Only two classes in the (K+1)-classes solution 
# are assigned the same class in the K-classes solution: the two which 
# are merged at this step 
output$combiM[[5]] %*% c(0,0,0,1,0,0) 
# recover the 5-classes soft clustering from the 6-classes soft clustering 
# and the 6 -> 5 combining matrix
all( output$combiz[[5]] == t( output$combiM[[5]] \%*\% t(output$combiz[[6]]) ) ) 
# the hard clustering under the 5-classes solution
head( output$classification[[5]] )
}
\keyword{ cluster }
