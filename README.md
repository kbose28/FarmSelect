
<!-- README.md is generated from README.Rmd. Please edit that file -->
FarmSelect: Factor Adjusted Robust Model Selection
==================================================

Goal of the package
-------------------

This R package implements a consistent model selection strategy for high dimensional sparse regression when the covariate dependence can be reduced through factor models. By separating the latent factors from idiosyncratic components, the problem is transformed from model selection with highly correlated covariates to that with weakly correlated variables. It is appropriate for cases where we have many variables compared to the number of samples. Moreover, it implements a robust procedure to estimate distribution parameters wherever possible, hence being suitable for cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature. See the paper on this method, Fan et al.(2017) <https://arxiv.org/pdf/1612.08490.pdf>, for detailed description of methods and further references.

Let there be *p* covariates and *n* samples. Let us model the relationship between the response vector *Y* and the covariates *X* as ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY%7D%20%3D%20%5Cmathbf%7BX%7D%5Cbm%7B%5Cbeta%7D%20+%20%5Cepsilon). Here *β* is a vector of size *p*. Non-zero values in this vector *β* denote which covariates truly belong in the model. For the covariates, assume the approximate factor model: ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX%7D%20%3D%20%5Cmathbf%7BFB%7D%5E%7BT%7D%20+%20%5Cmathbf%7BU%7D), where *F* are the *K* underlying factors, *B* are the factor loadings and *U* are the errors.

Installation
------------

You can install FarmSelect from github with:

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

There are two functions available.

-   `farm.select`: The main function farm.test which carries out the entire model testing procedure.
-   `farm.adjust`: Adjusts the data for latent fators.
-   `farm.mean`: Multivariate mean estimation with Huber's loss.
-   `farm.cov`: Multivariate covariance estimation with Huber's loss.

Also see the [`farm.scree`](https://www.rdocumentation.org/packages/FarmTest/versions/1.0.0/topics/farm.scree) function in the [`FarmTest`](https://cran.r-project.org/web/packages/FarmTest/index.html) package for how to generate diagnostic plots and output to evaluate the factor adjustment step.

Main function example: model selection
--------------------------------------

Here we generate data from a factor model with 3 factors. We have 50 samples of 100 dimensional data. The model is of size 5, where the first 5 covariates model coefficients being non-zero and the rest zero. The factors, loadings, erros are all generated from a normal distribution.

``` r
library(FarmSelect)
set.seed(100)
P = 100 #dimension
N = 50 #samples
K = 3 #nfactors
Q = 5 #model size
Lambda = matrix(rnorm(P*K, 0,1), P,K)
F = matrix(rnorm(N*K, 0,1), N,K)
UU = matrix(rnorm(P*N, 0,1), P,N)
X = Lambda%*%t(F)+UU
X = t(X)
beta_1 = 3+3*runif(Q)
beta = c(beta_1, rep(0,P-Q))
eps = rnorm(N)
Y = X%*%beta+eps 
output = farm.select(Y,X)
#> Call:
#> farm.select(Y = Y, X = X)
#> 
#>  Factor Adjusted Robust Model Selection 
#> loss function used: mcp
#> 
#> p = 100, n = 50
#> factors found: 3
#> size of model selected:
#>  5
```

``` r
names(output)
#> [1] "beta.chosen" "coef.chosen" "nfactors"    "X.res"       "Y.res"
```

The values X.res and Y.res are the covariates and responses after adjusting for latent factors. The formulas for these are ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY.res%7D%20%3D%20%28%5Cmathbf%7BI%7D_n-%5Cmathbf%7BP%7D%29%20%5Cmathbf%7BY%7D) and ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX.res%7D%20%3D%20%28%5Cmathbf%7BI%7D_n-%5Cmathbf%7BP%7D%29%20%5Cmathbf%7BX%7D%5E%7BT%7D), where ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BP%7D%3D%5Chat%7B%5Cmathbf%7BF%7D%7D%28%5Chat%7B%5Cmathbf%7BF%7D%7D%5E%7BT%7D%5Chat%7B%5Cmathbf%7BF%7D%7D%29%5E%7B-1%7D%5Chat%7B%5Cmathbf%7BF%7D%7D%5E%7BT%7D).

Now we use a different loss function for the model selection step.

``` r
output = farm.select(Y,X, loss = "lasso" )
#> Call:
#> farm.select(Y = Y, X = X, loss = "lasso")
#> 
#>  Factor Adjusted Robust Model Selection 
#> loss function used: lasso
#> 
#> p = 100, n = 50
#> factors found: 3
#> size of model selected:
#>  5
```

We may also use robust estimates of all the parameters. This may take more time, depending upon the dimension size of the problem.

``` r
output = farm.select(Y,X, robust = TRUE )
#> Call:
#> farm.select(Y = Y, X = X, robust = TRUE)
#> 
#>  Factor Adjusted Robust Model Selection 
#> loss function used: mcp
#> 
#> p = 100, n = 50
#> factors found: 3
#> size of model selected:
#>  5
```

Other functions
---------------

The function `farm.adjust` adjusts the dataset for latent factors. The number of factors is estimated internally by using the method in (Ahn and Horenstein 2013).

``` r
output = farm.adjust(Y,X)
names(output)
#> [1] "X.res"    "nfactors" "Y.res"
```

If known, we can provide this function (or the main function) the number of latent factors. Providing too large a number results in a warning message. The maximum number of factors possible is max(*n*, *p*) but a much smaller number is recommended.

``` r
output = farm.adjust(Y,X, K.factors  = 30)
#> Warning in farm.adjust(Y, X, K.factors = 30): Number of factors supplied is
#> > min(n,p)/2. May cause numerical inconsistencies
```

We see a warning telling us that it is not a good idea to calculate 30 eigenvalues from a dataset that has only 50 samples.

Notes
-----

1.  Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.

2.  The covariates do not need to be de-meaned before insertion into the function, this is done internally.

Ahn, SC, and AR Horenstein. 2013. “Eigenvalue Ratio Test for the Number of Factors.” Econometrica.
