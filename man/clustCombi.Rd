\name{clustCombi}
\alias{clustCombi}
\title{
Combining Gaussian Mixture Components for Clustering
}
\description{
Provides a hierarchy of combined clusterings from the EM/BIC Gaussian mixture solution to one class, following the methodology proposed in the article cited in the references.
}
\usage{
clustCombi(data, MclustOutput = NULL, ...)
}
\arguments{
  \item{data}{
		A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
	}
  \item{MclustOutput}{
		A list giving the optimal (according to BIC) parameters, conditional probabilities z, and loglikelihood, together with the associated classification and its uncertainty, as returned by Mclust. Please see Mclust documentation for the details of the components. Default value is \code{NULL}, in which case Mclust is called to get this output.
	}
  \item{\dots}{
		Optional arguments to be passed to called functions. Notably, any option (such as the numbers of components for which the BIC is computed; the models to be fitted by EM; initialization parameters for the EM algorithm...) to be passed to Mclust in case \code{MclustOutput} is \code{NULL}. Please see the \code{\link{Mclust}} documentation for more details.
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
\examples{

data(Baudry_etal_2010_JCGS_examples)
output <- clustCombi(ex4.1) # will run Mclust to get the MclustOutput

MclustOutput <- Mclust(ex4.1) # or you can run Mclust yourself
output <- clustCombi(ex4.1, MclustOutput) # and provide the output to clustCombi

# any further optional argument is passed to Mclust (see the Mclust documentation)
\dontrun{output <- clustCombi(ex4.1, modelName = "EEV", G = 1:15)}

output # is of class clustCombi
# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output, ex4.1) 

# the selected model and number of components by Mclust, ie by BIC with MLE 
# on Gaussian mixtures
output$MclustOutput 
# the selected number of components by Mclust: the combined hierarchy then 
# starts from this number of classes and ends at one
output$MclustOutput$G 
# the matrix whose [i,k]th entry is the probability that observation i in 
# the data belongs to the kth class according to the BIC solution
head( output$combiz[[output$MclustOutput$G]] ) 
# is the matrix whose [i,k]th entry is the probability that observation i in 
# the data belongs to the kth class according to the first combined 
# ((output$MclustOutput$G-1)-classes) solution
head( output$combiz[[output$MclustOutput$G-1]] ) 
# the matrix describing how to merge the 6-classes solution to get the 
# 5-classes solution
output$combiM[[5]] 
# for example the following code returns the label of the class (in the 
# 5-classes combined solution) to which the 4th class (in the 6-classes
# solution) is assigned. Only two classes in the (K+1)-classes solution 
# are assigned the same class in the K-classes solution: the two which 
# are merged at this step... 
output$combiM[[5]] %*% c(0,0,0,1,0,0) 
# recover the 5-classes soft clustering from the 6-classes soft clustering 
# and the 6 -> 5 combining matrix
all( output$combiz[[5]] == t( output$combiM[[5]] \%*\% t(output$combiz[[6]]) ) ) 
# the hard clustering under the 5-classes solution
head( output$classification[[5]] )
}
\keyword{ cluster }
