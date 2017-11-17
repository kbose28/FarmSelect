
<!-- README.md is generated from README.Rmd. Please edit that file -->
FarmSelect: Factor Adjusted Robust Model Selection
==================================================

Goal of the package
-------------------

This R package implements a consistent model selection strategy for high dimensional sparse regression when the covariate dependence can be reduced through factor models. By separating the latent factors from idiosyncratic components, the problem is transformed from model selection with highly correlated covariates to that with weakly correlated variables. It is appropriate for cases where we have many variables compared to the number of samples. Moreover, it implements a robust procedure to estimate distribution parameters wherever possible, hence being suitable for cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature. See the paper on this method, Fan et al.(2017)<https://arxiv.org/pdf/1612.08490.pdf>, for detailed description of methods and further references.

The observed data *x*<sub>*i*, *j*</sub> is assumed to follow a factor model ![equation](https://latex.codecogs.com/gif.latex?\mathbf%7Bx%7D_i&space;=&space;\mathbf%7BBf%7D_i&space;+\mathbf%7Bu%7D_i), where *f* are the underlying factors, *B* are the factors loadings, *u* are the errors, and *μ* is the mean effect to be tested. We assume the data is of dimension *p* and the sample size is *n*, leading to *p* hypothesis tests.

Installation
------------

You can install FarmTest from github with:

``` r
install.packages("devtools")
devtools::install_github("kbose28/FarmSelect")
library(FarmSelect)
```

Getting help
------------

Help on the functions can be accessed by typing "?", followed by function name at the R command prompt.

Issues
------

-   Error: "...could not find build tools necessary to build FarmSelect": Since `FarmSelect` relies on `C++` code, command line tools need to be installed to compile the code. For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See (<https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites>).

-   Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue.

    1.  In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions here: <http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/>.

    2.  For &gt;= R 3.4.\* : download the installer from the here: <https://gcc.gnu.org/wiki/GFortranBinaries#MacOS>. Now simply run the installer. (If installer is not available for your version of OS, use the latest one.)

-   Error: "... .rdb': No such file or directory" Try devtools::install\_github("kbose28/FarmSelect", dependencies=TRUE)

-   Error in RStudio even after installing XCode: "Could not find tools necessary to build FarmSelect": This is a known bug in RStudio. Try options(buildtools.check=function(action) TRUE) in RStudio to prevent RStudio from validating build tools.

Functions
---------

There are three functions available.

-   `farm.test`: The main function farm.test which carries out the entire hypothesis testing procedure.
-   `farm.FDR`: Apply FDR control to a list of input p-values. This function rejects hypotheses based on a modified Benjamini- Hochberg procedure, where the proportion of true nulls is estimated using the method in \[@storey2015\].
-   `farm.scree`: Estimate the number of factors if it is unknown. The farm.scree function also generates two plots to illustrate how the number of latent factors is calculated.

Simple hypothesis testing example
---------------------------------

Here we generate data from a factor model with 3 factors. We have 20 samples of 100 dimensional data. The first five means are set to 2, while the other ones are 0. We conduct a hypotheses test for these means.

``` r
#library(FarmSelect)
#set.seed(100)
#p = 100
#n = 20
#epsilon = matrix(rnorm( p*n, 0,1), nrow = n)
#B = matrix(rnorm(p*3,0,1), nrow=p)
#fx = matrix(rnorm(3*n, 0,1), nrow = n)
#mu = rep(0, p)
#mu[1:5] = 2
#X = rep(1,n)%*%t(mu)+fx%*%t(B)+ epsilon
#output = farm.test(X)
```

Now we carry out a one-sided test, with the FDR to be controlled at 1%. Then we examine the output

``` r
#output = farm.test(X, alpha = 0.01,alternative = "greater")
#names(output)
#print(output$rejected)
#hist(output$means, 10, main = "Estimated Means", xlab = "")
```

Other functions
---------------

The function `farm.scree` makes some informative plots. It is possible to specify the maximum number of factors to be considered and the maximum number of eigenvalues to be calculated in this function. We recommend min(n,p)/2 as a conservative threshold for the number of factors; this also prevents numerical inconsistencies like extremely small eigenvalues which can blow up the eigenvalue ratio test.

``` r
#output = farm.scree(X, K.factors = 15, K.scree = 10)
```

We see a warning telling us that it is not a good idea to calculate 15 eigenvalues from a dataset that has only 20 samples.

Let us generate data from a Gaussian distribution with mean 0. Suppose we perform a simple `t.test` in R and need to adjust the output p-values for multiple testing. The function `farm.FDR` lets us carry out multiple comparison adjustment and outputs rejected hypotheses. We see that there are no rejections, as expected from a zero-mean Gaussian distribution.

``` r
#set.seed(100)
#Y = matrix(rnorm(1000, 0, 1),10)
#pval = apply(Y, 1, function(x) t.test(x)$p.value)
#output = farm.FDR(pval)
#output$rejected
```

Notes
-----

1.  If some of the underlying factors are known but it is suspected that there are more confounding factors that are unobserved: Suppose we have data ![equation](https://latex.codecogs.com/gif.latex?X%20%3D%20%5Cmu%20+%20Bf%20+%20Cg%20+%20u), where *f* is observed and *g* is unobserved. In the first step, the user passes the data {*X*, *f*} into the main function. From the output, let us construct the residuals: ![equation](https://latex.codecogs.com/gif.latex?Xres%20%3D%20X%20-%20Bf). Now pass ![equation](https://latex.codecogs.com/gif.latex?Xres) into the main function, without any factors. The output in this step is the final answer to the testing problem.

2.  Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.

3.  The farm.FDR function uses code from the [`pi0est`](https://www.rdocumentation.org/packages/qvalue/versions/2.4.2/topics/pi0est) function in the [`qvalue`](http://bioconductor.org/packages/release/bioc/html/qvalue.html) package \[@storey2015\] to estimate the number of true null hypotheses, and inherits all the options from `pi0est`.

4.  See individual function documentation for detailed description of methods and their references.
