
MCLUST is a software package for cluster analysis based on parameterized 
Gaussian mixture models.

MCLUST is written in Fortran with an interface to the S-PLUS commercial 
package (http://www.mathsoft.com/splus). Commercial S-PLUS includes two
functions `mclust' and `mclass', which have been superseded by functions 
`mhtree' and `mhclass' in MCLUST. In addition, MCLUST a great many other
functions for cluster analysis (see below).

copyright 1996, 1998 Department of Statistics, University of Washington
funded by ONR contracts N00014-96-1-0192 and N00014-96-1-0330

Permission granted for unlimited redistribution for non-commercial use only.

Please report all anomalies to fraley@stat.washington.edu.

References : (see http://www.stat.washington.edu/fraley/reports.shtml)

   Banfield and Raftery, Biometrics 49:803-821, 1993.

   DasGupta and Raftery, J. Amer. Stat. Assoc. 93:294-302, 1998.
      (earlier version Tech. Rep. 295)

   Fraley, SIAM J. on Scientific Computing 20(1):270-281, 1999.
      (earlier version Tech. Rep. 311)

   Fraley and Raftery, Tech. Rep. 329, U. of Washington Statistics, Feb. 1998.
      (to appear in The Computer Journal, 
       special issue on Clustering and Classification)

   Fraley and Raftery, Tech. Rep. 342, U. of Washington Statistics, Nov. 1998.

-----------------------------------------------------------------------------

Anyone can use MCLUST for noncommercial purposes, but in order to do so 
without modification they must have access to the S-PLUS (version 3.3 or
higher).

The most up to date version of MCLUST can be obtained via anonymous ftp from 
ftp.u.washington.edu in the directory public/mclust. There is a gzipped tar
file `MCLUST.tar.gz' for UNIX users, and a WinZip file `mclust.zip' for
Windows users.

---------------------------------------------------------------------------
MCLUST Installation:
---------------------------------------------------------------------------

UNIX:
-----

The file MCLUST.tar.gz is a compressed version of a directory called MCLUST
containing the MCLUST software and documentation. The commands for extracting
this directory (and removing the tar file) are:

%gunzip MCLUST.tar.gz
%tar xvf MCLUST.tar
%rm MCLUST.tar  

The following commands compile the Fortran code to create an object file
for dynamic or static loading into S-PLUS:

%cd MCLUST
%Splus COMPILE mclust.o

To statically load the object file into S-PLUS, move it to the working 
directory and use the following command:

%Splus LOAD mclust.o

This creates a file called `local.Sqpe', a local version of S-PLUS that 
includes the MCLUST Fortran functions. When S-PLUS is invoked in a directory
containing `local.Sqpe', this local version of S-PLUS will be run. 
Dynamic loading, which is available in S-PLUS on some but not all UNIX 
platforms, may be preferable to static loading since `local.Sqpe' is a
large file. 

The following commands put the S-PLUS functions in the MCLUST/.Data directory:

%cd MCLUST
%cat mclust.S | Splus

Help files are located in MCLUST/.Data/.Help. You can access eveything from 
your S-PLUS session if you make MCLUST a subdirectory of your working 
directory and invoke the command

> attach("MCLUST/.Data")

in S-PLUS. You can put this `attach' command in a `.First' function if you 
want the MCLUST S_PLUS functions to be attached automatically when you enter 
S-PLUS.

Windows:  
--------

The file `mclust.zip' is a compressed version of a folder called `mclust' 
containing the MCLUST documentation, S-PLUS functions, help files, and 
WATCOM Fortran object code. Open `mclust.zip' to invoke WinZip, and then
extract the `mclust' folder. In what follows we assume the `mclust' folder
has been placed in your working Splus directory.

You can put the S-PLUS functions for MCLUST into your working `_Data'
directory with the S-PLUS command:

> source("mclust/mclust.S")

Alternatively, you may want to place them in a separate directory to and use
the `attach' command in S-PLUS to access them, to ensure that you don't
change them by accident.

Dynamically load the object code `mclust.obj' with the S-PLUS command:

> dyn.load("mclust/mclust.obj")

This command  needs to be executed each time you enter S-PLUS. You may
want to include this command in a `.First' file if you plan to use MCLUST 
every time you run S-PLUS in this directory.

The help files are located in the `mclust.HLP' archive.

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

> iris.matrix <- matrix(aperm(iris,c(1,3,2)),150,4,dimn=dimnames(iris)[1:2])

The function mhtree() for hierarchical clustering can be applied to this data 
to obtain a classification tree:

> cltree <- mhtree(iris.matrix)

The classification produced by mhtree() for various numbers of clusters
can be obtained via the function mhclass(). For example, to obtain the
classification for 2 through 4 clusters:

> cl <- mhclass(cltree, 2:4)

which produces a matrix whose columns represent the mhtree() classification
for 2, 3, and 4 clusters, respectively.

A given classification can be plotted using the function clpairs():

> clpairs(iris.matrix, cl[,"2"])
> clpairs(iris.matrix, cl[,"3"])

plots the 3-cluster class from mhtree().
The number of dimensions may need to be reduced in order to get a good
display of the classification.

mhtree() will accept an initial coarse partition as input if one is
available.

In the above example, the default model (Gaussian mixtures with no
constraints on the means and covariance matrices of each component).
Other parameterizations of Gaussian mixtures are available in mhtree().
They differ as to whether certain geometric characteristics (volume, shape, 
and or orientation) are uniform among clusters.

For example, the following produces a classification tree for hierarchical 
clustering based on a spherical covariances with equal volumes

> cltree <- mhtree( iris.matrix, modelid = "EI")
> clpairs(iris.matrix, mhclass(cltree,2))

This is equivalent to the sum-of-squares method or Ward's method.

Other model options that are available for mhtree() are described in the
help file.

Classifications produced by mhtree() are often good but not necessairly
optimal. MCLUST provides iterative `EM' (Expectation-Maximization) methods 
for maximum likelihood clustering with parameterized Gaussian mixture models. 
`EM' iterates between an `E'-step, which computes a matrix `z' whose ik-th
element is an estimate of the conditional probability that observation i
belongs to group k given the current parameter estimates, and an `M-step', 
which computes maximum-likelihood parameter estimates given $z$. 
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
For example, me() can be used to `refine' a classification produced by
mhtree():

> cltree <- mhtree( iris.matrix, modelid = "VI") # non-uniform spherical model
> cl <- mhclass(cltree, 3) # 3-group mhtree classification
> z <- me( iris.matrix, modelid = "VI", ctoz(cl)) # optimal z

The function ctoz() converts a discrete classification into
a `z' matrix that has only 0,1 entries with exactly one 1 per row.
In general, the models used in mhtree() and me() need not be the same.
It may be desirable to use one of the faster methods in mhtree()
(e.~g. spherical or unconstrained models), followed by specification of
a more complex model in EM.

To plot the classification produced by me(), do

> clpairs( iris.matrix, ztoc(z))

The function ztoc() converts conditional probabilties into a discrete
classification by classifying each observation into the class corresponding
to the column in which the `z' value for that observation is maximized.

The uncertainty in the classification corresponding to `z' can be obtained
as follows:

> uncer <- 1 - apply( z, 1, "max")

The Splus function quantile() applied to uncertainly gives a measure of the 
quality of the classification. 

> quantile(uncer)

In this case the indication is that the majority of observations are well
classified. Note, however, that when groups intersect
uncertain classifications would be expected in the overlapping regions.

The function mstep() allows recovery of the parameters in the model:

> pars <- mstep( iris.matrix, modelid = "VI", z)
> pars

The means and standard deviation of the clusters (and optionally the
partition or classification) may be plotted using the function mixproj(),
which 

> mixproj(iris.matrix, ms=pars, part=ztoc(z), dimens = c(1,2))
> mixproj(iris.matrix, ms=pars, part=ztoc(z), dimens = c(3,4))

It's also possible to compute the Bayesian Information Criterion, or `BIC'
for a given data set, model, and conditional probability estimates.

> cltree <- mhtree(iris.matrix)  # uses default unconstrained model
> cl <- mhclass(cltree, 2:3)
> z <- me( iris.matrix, modelid = "EI", ctoz(cl[,"2"])) # uniform spherical
> bic(iris.matrix, modelid = "EI", z)
> z <- me( iris.matrix, modelid = "EI", ctoz(cl[,"3"]))
> bic(iris.matrix, modelid = "EI", z)
> z <- me( iris.matrix, modelid = "VI", ctoz(cl[,"2"])) # spherical
> bic(iris.matrix, modelid = "VI", z)
> z <- me( iris.matrix, modelid = "VI", ctoz(cl[,"3"]))
> bic(iris.matrix, modelid = "VI", z)

The function unclass() applied to the bic value will reveal an attribute 
labeled {\tt "rcond"} representing the reciprocal condition estimate will, 
which falls between 0 and 1. 
The closer this estimate is  to 0, the more likely it is that the estimated 
values are numerically inaccurate due to a covariance that is nearly singular.
You should check this quantity if you are getting bic values that don't
seem to make sense.

Of these four possiblities (2 models x 2 numbers of groups) illustrated above,
the best model and number of clusters is "VI" (spherical), with 3 groups.

MCLUST provides two functions, emclust1() and emclust(), for cluster analysis 
with BIC. The input to emclust1() is the data, desired numbers of groups,
and a pair of models, one to be used in the hierarchical clustering phase,
and the other to be used in the EM phase:

> bicvals <- emclust1( iris.matrix, nclus = 1:6, modelid = c("VI", "VVV"))
> bicvals
> plot(bicvals)
> sumry <- summary(bicvals, iris.matrix) # summary object for emclust1()
> sumry

The input to emclust() is data, desired numbers of groups, and a list of
models to apply in the EM phase (it uses the unconstrained model for
hierarchical clustering). By default, all available models are used:

> bicvals <- emclust( iris.matrix, nclus = 1:6)
> bicvals
> plot(bicvals)
> sumry <- summary(bicvals, iris.matrix) # summary object for emlust()
> sumry

The best model among those fitted by emclust() is the uniform shape model 
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

The EM functions in MCLUST have the option of including a Poisson term for 
noise in the model. For me(), mstep(), estep() and bic(), use the `noise = T' 
option to use this model. The initial conditional probabilities for the noise 
term correspond to the last column of the conditional probabilities.

In emclust1() and emclust(), an initial estimate of the noise (in the form of
a logical vector indicating which observations are to be considered noise) 
must be supplied in order to specify this model. The Splus function NNclean() 
in statlib is one option for estimating noise --- see Byers and Raftery, 
tech rep no 305, U of Washington Statistics (April 1996) --- to appear in JASA
in June 1998). To obtain NNclean, send a message of the form
`send nnclean from S' to statlib@lib.stat.cmu.edu.

The following reproduces the `chevron' minefield example from
Fraley and Raftery (1998) for which the data is included in this distribution:

> chevron.nnclean <- NNclean(chevron, k = 5) # 5 nearest-neighbor denoising
> chevron.noise <- !chevron.nnclean$z
> chev.bic <- emclust(chevron, noise = chevron.noise)

------------------------------------------------------------------------------
