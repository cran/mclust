                        MCLUST/EMCLUST Software
                                 by 
                    Chris Fraley and Adrian Raftery 
              Department of Statistics, University of Washington

URL:             http://www.stat.washington.edu/mclust

Sponsors:

                  Office of Naval Research 
        (contracts N00014-96-1-0192 and N00014-96-1-0330)
                 
                 National Institutes of Health 
                 (contract 1 R01 CA94212-01)

MCLUST is a software package for clustering, discriminant analysis and
density estimation. It includes hierarchical clustering and EM for 
parameterized Gaussian mixture models, and has functions that incorporate 
these into a comprehensive clustering strategy which use the Bayesian 
Information Criterion (BIC) to select the model and number of clusters.
It also includes functions for discriminant analysis and density estimation.

MCLUST is written in Fortran with an interface to the S-PLUS commercial 
package (http://www.splus.mathsoft.com). Commercial S-PLUS includes two
functions `mclust' and `mclass', which have been superseded by functions 
`hc' and `hclass' in MCLUST. 

copyright 1996, 1998, 2002 Department of Statistics, University of Washington
funded by ONR contracts N00014-96-1-0192 and N00014-96-1-0330
and NIH Grant 1 R01 CA94212-01.

Permission granted for unlimited redistribution for non-commercial use only.
Commerical use requires a licensing agreement with the University of 
Washington.

Please report all anomalies to fraley@stat.washington.edu.

References : (see http://www.stat.washington.edu/mclust)

   C. Fraley and A. E. Raftery, MCLUST: Software for model-based clustering, 
   discriminant analysis, and density estimation, 
   Technical Report No. 415, Department of Statistics, 
   University of Washington, September 2002.

   C. Fraley and A. E. Raftery, Model-based clustering, discriminant analysis,
   and density estimation, Journal of the American Statistical Association
   97:611-631 (2002)

   C. Fraley, Algorithms for model-based Gaussian hierarchical
   clustering, SIAM Journal on Scientific Computing 20:270-281 (1998).

-----------------------------------------------------------------------------

Anyone can use MCLUST for noncommercial purposes, but in order to do so 
without modification they must have access to the S-PLUS commercial
package (version 6.0 or higher for UNIX or Windows). [Do not use
on Compaq Alpha]

The most up-to-date version of MCLUST can be obtained via anonymous ftp from 
ftp.u.washington.edu in the directory public/mclust. There is a gzipped tar
file `MCLUST.tar.gz' for UNIX users, and a WinZip file `mclust.zip' for
Windows users.

---------------------------------------------------------------------------
MCLUST Installation:
---------------------------------------------------------------------------

UNIX (S-PLUS ):
---------------

(Note: MCLUST does not work on the Compaq Alpha platform)

The file MCLUST.tar.gz is a compressed version of a directory called MCLUST
containing the MCLUST software and documentation. The commands for extracting
this directory (and removing the tar file) are:

%gunzip MCLUST.tar.gz
%tar xvf MCLUST.tar
%rm MCLUST.tar  

The following commands compile the Fortran code to create an object file
for dynamic or static loading into S-PLUS:

%cd MCLUST
%Splus CHAPTER
%Splus make
%Splus TRUNC_AUDIT 0

The `Splus make' command 
- compiles the Fortran code, which will automatically loaded when S-PLUS is 
invoked in this directory.
- sources the S-PLUS functions in file mclust.q into the {\tt .Data} directory
- compiles and loads the help (.sgml) files 
The software can now be used by running S-PLUS in directory MCLUST. 
They can also be run in a different working directory,
using the following S-PLUS commands
> dyn.open(".../MCLUST/S.so")  
> attach(".../MCLUST/.Data")    
where .../ indicates the path to the directory MCLUST.

Windows (S-PLUS 6):
----------------------------

The file `mclust.zip' contains the necessary files for running MCLUST
with S-PLUS 6 for Windows. To install MCLUST on Windows, unzip this file
to obtain a directory called `mclust'. That directory needs to be 
attached in S-PLUS in order to access the MCLUST functions. 
To use the MCLUST help files, either click on the file `mclust.chm'
included in the folder, or else execute it as a command from a DOS prompt.

---------------------------------------------------------------------------
MCLUST Usage:
---------------------------------------------------------------------------

To illustrate the use of MCLUST  we use the `iris' data in S-PLUS for 
which a description is available from the command line via `help(iris)' or 
from the help window under `datasets'. The iris data is a 3-dimensional array,
giving 4 measurements on 50 flowers from each of 3 species of iris. 
MCLUST  requires 2-dimensional numeric data as input, so we collapse
the iris data into a matrix in which row dimension corresponds to the 150
individual flowers and the column dimension corresponds to the 4 measurements.
The species information has been discarded.

> irisMatrix <- matrix(aperm(iris,c(1,3,2)),150,4,dimn=dimnames(iris)[1:2])

The function hc() for hierarchical clustering can be applied to this data 
to obtain a classification tree:

> cltree <- hc(irisMatrix)

The classification produced by hc() for various numbers of clusters
can be obtained via the function hclass(). For example, to obtain the
classification for 2 through 4 clusters:

> cl <- hclass(cltree, 2:4)

which produces a matrix whose columns represent the hc() classification
for 2, 3, and 4 clusters, respectively.

A given classification can be plotted using the function clPairs():

> clPairs(irisMatrix, cl[,"2"])
> clPairs(irisMatrix, cl[,"3"])

plots the 3-cluster class from hc().
The number of dimensions may need to be reduced in order to get a good
display of the classification.

hc() will accept an initial partition as input if one is available.

In the above example, the default model (Gaussian mixtures with no
constraints on the means and covariance matrices of each component).
Other parameterizations of Gaussian mixtures are available in hc().
They differ as to whether certain geometric characteristics (volume, shape, 
and or orientation) are uniform among clusters.

For example, the following produces a classification tree for hierarchical 
clustering based on a spherical covariances with equal volumes

> cltree <- hc( irisMatrix, modelName = "EII")
> clPairs(irisMatrix, hclass(cltree,2))

This is equivalent to the Ward's sum-of-squares  method.

MCLUST also provides iterative `EM' (Expectation-Maximization) methods 
for maximum-likelihood clustering with parameterized Gaussian mixture models. 
`EM' iterates between an `E'-step, which computes a matrix `z' whose ik-th
element is an estimate of the conditional probability that observation i
belongs to group k given the current parameter estimates, and an `M-step', 
which computes maximum-likelihood parameter estimates given `z'. 
Basically, `z' can be thought of as a continuous representation of a 
classification; `z' has a row for each observation and an column for each 
group, and its entries entries of `z' are lie in the interval [0,1],
with each row summing to 1. 
In limit, the parameters converge to the maximum likelihood values
for the Gaussian mixture model, and the sums of the columns of the
`z' converge to `n' times the mixing proportions, where `n' is the
number of observations in the data.

MCLUST provides a function me() that does the M-step followed by the
E-step for EM. Given the data, initial estimates for the conditional 
probabilities `z',  and the model specification, it produces the conditional 
probabities associated with the maximum likelihood parameters.
Initial estimates of `z' may be obtained from a discrete classification.
For example, me() can be used to `refine' an available classification 
such as one produced by hc():

> cltree <- hc( irisMatrix, modelName = "VII") # non-uniform spherical model
> cl <- hclass(cltree, 3) # 3-group hc classification
> irisModel <- me( irisMatrix, modelName = "VII", unmap(cl)) 

The function unmap() converts a discrete classification into
a `z' matrix that has only 0,1 entries with exactly one 1 per row.
In general, the models used in hc() and me() need not be the same.
It may be desirable to use one of the faster methods in hc()
(e.~g. spherical or unconstrained models), followed by specification of
a more complex model in EM.

To plot the classification produced by me(), do

> z <- irisModel$z ## optimal z produced my me()
> clPairs( irisMatrix, map(z))

The function map() converts conditional probabilties into a discrete
classification by classifying each observation into the class corresponding
to the column in which the `z' value for that observation is maximized.

To plot coordinate or random projections, do

> do.call("coordProj", c(list(data = irisMatrix, dimens = c(1,2), scale = T), 
          irisModel))
> do.call("randProj", c(list(data = irisMatrix, seed = 0, scale = T), 
          irisModel))

Note that if you don't use the scale = T option, it won't be obvious that
a spherical model is being fitted.

The uncertainty in the classification corresponding to `z' can be obtained
as follows:

> uncer <- 1 - apply( z, 1, "max")

The Splus function quantile() applied to uncertainly gives a measure of the 
quality of the classification. 

> quantile(uncer)

In this case the indication is that the majority of observations are well
classified. Note, however, that when groups intersect
uncertain classifications would be expected in the overlapping regions.

Parameters in the model and other infirmation can also be recovered from
the me() output:

> names(irisModel)
 [1] "n"         "d"         "G"         "z"         "mu"        "sigma"    
 [7] "sigmasq"   "decomp"    "pro"       "loglik"    "Vinv"      "modelName"

MCLUST also provides a function em() for iterated E-step followed by
M-step, as well as estep() and mstep() for the individual steps of
EM, for each of the models.

It is also possible to compute the Bayesian Information Criterion, or `BIC'
for a given data set, and optimal model parameters. For the above example,
this is accompliched by:

> do.call( "bic", irisModel)

MCLUST provides a function EMclust() for cluster analysis with BIC. 
The input to EMclust() is data, desired numbers of groups, a list of
models to apply in the EM phase, and a set of merge pairs from
hierarchical clustering (it computes these for model-based hierarchical
clustering via the unconstrained model by default). 
All available models are used with up to 9 groups are tried unless otherwise
specified.

> irisBIC <- EMclust(irisMatrix)
> irisBIC
> plot(irisBIC)
> sumry <- summary(irisBIC, irisMatrix) # summary object for EMclust()
> sumry

The best model among those fitted by EMclust() is the uniform shape model 
"VEV", with 2 clusters. The same model with 3 clusters has a BIC value that is
not signifcantly different from the maximum, so that the conclusion is that
there are either 2 or 3 clusters in the data. The 2 cluster EM result 
separates the first species from the other two, while the 3 cluster result 
nearly separates the three species (there are 5 missclassifications out of 
150). 

The summary function allows the summary to be taken over a subset of the
number of clusters (as well as models in the case of emclust()) for which 
the analysis was originally run.
It returns a list whose components are the parameter and z values
for the best fit among the models selected according to BIC, as well
as the associated classification and its uncertainty.

For a complete analysis, it may be desirable to try different permutations of 
the observations, use different subsets of the observations, and/or perturb 
the data, to see if the classification remains stable.
Scaling or otherwise transforming the data may also affect the results.
One should always examine the data beforehand, in case the dimensions can
be reduced due to highly correlated variables.

------------------------------------------------------------------------------
POISSON NOISE
------------------------------------------------------------------------------

The EM functions in MCLUST have the option of including a Poisson term for 
noise in the model. For me(), em(), mstep(), estep() and bic(), use the 
`noise = T' option to use this model. The initial conditional probabilities 
for the noise term correspond to the last column of the conditional 
probabilities.

The function EMclustN() is provided for model-based clustering with Poisson
noise. An initial estimate of the noise (in the form of a logical vector 
indicating which observations are to be considered noise) must be supplied in 
order to specify this model. The S-PLUS function NNclean() 
in StatLib (www.lib.stat.cmu.edu) is one option for estimating noise --- 
see Byers and Raftery, Nearest-neighbor clutter removal for estimating 
features in spatial point processes, Journal of the American Statistical
Association, 93:577-584, June 1998 (earlier version Tech. Rep.
no. 305). To obtain NNclean, send a message of the form `send nnclean from S' 
to  statlib@lib.stat.cmu.edu.

The following reproduces the `chevron' minefield example from Fraley and 
Raftery, How many clusters? Which clustering method? --- Answers via
Model-Based Clter Analysis,  Computer Journal 41:578-588 (1998) 
for which the data is included in this distribution:

> chevron.nnclean <- NNclean(chevron[,-1], k = 5) #5 nearest-neighbor denoising
> chevBIC <- EMclustN(chevron[,-1], noise = !chevron.nnclean$z)
> summary(chevBic)

The first 350 observations are in the minefield; the rest are noise.

------------------------------------------------------------------------------
DISCRIMINANT ANALYSIS
------------------------------------------------------------------------------

The function estep() in MCLUST can be used for discriminant analysis with the 
various parameterized Gaussian models. 
Below we illustrate some possible discriminant analysis techniques
using the Lansing Woods Data (D. J. Gerrard, Competition quotient: 
a new measure of the competition affecting individual forest trees. 
Researach Bulletin No. 20, Agricultural Experimental Station, Michigan State Univeristy, 1969) in which the object is to distinguish the hickory and 
maple trees. We use the odd-numbered observations as training data and 
the even-numbered observations as test data.
The first method uses `mstep' and `estep' for discriminant analysis,
and seletcs the model via BIC rather than cross-validation:

> lansing
> n <- nrow(lansing)

> odd <- seq(1, n, by = 2)

> oddBic <- rep(0, length(.Mclust$emModelNames))
> names(oddBic) <- .Mclust$emModelNames
> for (M in .Mclust$emModelNames) {
     oddMstep <- mstep(modelName = M, data = lansing[odd,-3], 
                                 z = unmap(lansing[odd,3]))
     oddEstep <- do.call("estep", c(list(modelName=M, data=lansing[odd,-3]), 
                                    oddMstep))
     oddBic[M] <- do.call("bic", oddEstep)
} 

> oddBic
> oddModelName <- (names(oddBic)[oddBic == max(oddBic)])[1]
> oddModel <- mstep(modelName = oddModelName, data = lansing[odd,-3], 
                                 z = unmap(lansing[odd,3]))

> oddEstep <- do.call("estep", c(list(data = lansing[,-3]), oddModel))

The following summarizes the clustering information:

data.frame(class = lansing[,3], oddEstep$z)

Alternatively, the density of the data can computed relative to a model fitted to each class:

> oddH <- lansing[odd,3] == "hickory"
> oddM <- lansing[odd,3] == "maple"
> hBic <- mBic <- rep(0, length(.Mclust$emModelNames))
> names(hBic) <- names(mBic) <- .Mclust$emModelNames
> for (M in .Mclust$emModelNames) {
     hMstep <- mstep(modelName = M, data = lansing[oddH,-3], 
                                 z = unmap(lansing[oddH,3]))
     mMstep <- mstep(modelName = M, data = lansing[oddM,-3], 
                                 z = unmap(lansing[oddM,3]))
     hEstep <- do.call("estep", c(list(modelName=M, data=lansing[oddH,-3]), 
                                       hMstep))
     mEstep <- do.call("estep", c(list(modelName=M, data=lansing[oddM,-3]), 
                                       mMstep))
     hBic[M] <- do.call("bic", hEstep)
     mBic[M] <- do.call("bic", mEstep)
} 

> hModelName <- (names(hBic)[hBic == max(hBic)])[1]
> mModelName <- (names(mBic)[mBic == max(mBic)])[1]
> hModel <- mstep(modelName = hModelName, data = lansing[oddH,-3], 
                                 z = unmap(lansing[oddH,3]))
> mModel <- mstep(modelName = mModelName, data = lansing[oddM,-3], 
                                 z = unmap(lansing[oddM,3]))

> hDens <- do.call("dens", c(list(data = lansing[,-3]), hModel))  
> mDens <- do.call("dens", c(list(data = lansing[,-3]), mModel))  

> data.frame(class = lansing[,3], h = hDens, m = mDens)

MCLUST also provides functions implementing a more powerful (but more
computationally intensive) discriminant analysis method called MclustDA
(Fraley and Raftery, Model-based clustering, discriminant analysis,
 and density estimation, Journal of the American Statistical Association
97:611-631, June 2002):

> train <- mclustDAtrain(data = lansing[odd,-3], labels = lansing[odd,3])
> test <- mclustDAtest(data = lansing[,-3], models = train)

> data.frame(class = lansing[,3], test)

-------------------------------------------------------------------------------
LARGE DATA SETS
-------------------------------------------------------------------------------

The discrimination capability of MCLUST can be used for classification of
large data sets. First, cluster analysis with the methodology of EMclust()
can be performed on a subset of the data. Then the remaining data points can 
then be classified (in reasonable sized blocks) using estep() or dens(), 
using the optimal parameters from the subset  as illustrated above.
EMclust() also includes a provision for using a subset of the data in the 
hierarchical clustering phase before applying EM to the full data set. 
This strategy is often adequate for large data sets.

-------------------------------------------------------------------------------
HIGH DIMENSIONAL DATA
-------------------------------------------------------------------------------

Models in which the orientation is allowed to vary between clusters (EEV, VEV,
and VVV in the current version of MCLUST) have O(d^2) parameters per cluster,
where d is the dimension of the data. For this reason, MCLUST may not work 
well or may otherwise be inefficient for these models when applied to 
high-dimensional data. It may still be possible to analyze such data with
MCLUST by restriction to models with fewer parameters (the spherical models 
EII, VII, diagonal models EEI, VEI, EVI, VVI,  and the constant variance model
EEE in MCLUST), or else by applying a dimension-reduction technique such as 
principal components. 

For datasets in which the number of observations is smaller than the data 
dimension, it may be possible to use the spherical or diagonal models.

