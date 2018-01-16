#' @useDynLib FarmSelect
#' @importFrom  Rcpp sourceCpp
#' @importFrom  graphics mtext plot points axis par barplot
#' @import methods
#' @import utils
#' @import grDevices
#' @import stats
#' @import glmnet
#' @import ncvreg
NULL
###################################################################################
## Main function that performs model selection
###################################################################################
#' Factor-adjusted robust model selection
#'
#' Given a covariate matrix and output vector, this function first adjusts the covariates for underlying factors and then performs model selection.
#' @param Y a size n outcome vector.
#' @param X an n x p covariate matrix with each row being a sample. Must have same number of rows as the size of \code{Y}.
#' @param loss an \emph{optional} character string specifying the loss function to be minimized. Must be one of "mcp" (default), "scad" or "lasso". You can just specify the initial letter.
#' @param verbose an \emph{optional} boolean determining whether or not to print output to the console. Default is TRUE.
#' @param robust an \emph{optional} boolean, specifying whether or not to use robust estimators for mean and variance. Default is FALSE.
#' @param \dots Arguments passed to the \code{\link{farm.adjust}} function.
#' @return A list with the following items
#' \item{beta.chosen}{the indices of the covariates chosen in the model}
#' \item{coef.chosen}{the coefficients of the chosen covariates}
#' \item{X.res}{the decorrelated covariate matrix}
#' \item{Y.res}{the response vector corresponding to the decorrelated covariates}
#' \item{nfactors}{number of (estimated) factors}
#' @details Number of rows and columns of the covariate matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For details about the method, see Fan et al.(2017).
#' @details For formula of how the covariates and the repsonse vector are  adjusted for latent factors, see Section 3.2 in Fan et al.(2017).
#' @examples
#' set.seed(100)
#' P = 100 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 5 #model size
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' beta_1 = 3+3*runif(Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps = rnorm(N)
#' Y = X%*%beta+eps
#' output = farm.select(Y,X)
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/abs/1612.08490}
#' @export
farm.select <- function (Y, X,  loss=c("mcp", "scad", "lasso"), verbose = TRUE, robust = FALSE, ...){
  #error checking
  X = t(X)
  if(NCOL(X)!=NROW(Y)) stop('number of rows in covariate matrix should be size of the outcome vector')
    output.final = farm.select.adjusted(Y, X,   loss=c("mcp", "scad", "lasso"), ...)
    if(verbose){output.final.call = match.call()
    cat("Call:\n")
    print(output.final.call)
    cat("\n Factor Adjusted Robust Model Selection \n")
    cat(paste("loss function used: ",  match.arg(loss), "\n",sep = ""))
    cat(paste("\np = ", NROW(X),", n = ", NCOL(X),"\n", sep = ""))
    cat(paste("factors found: ", output.final$nfactors, "\n", sep = ""))
    cat("size of model selected:\n")
    if(is.character(output.final$beta.chosen)){ cat(" no variable selected\n")} else{ cat(paste(" ", NROW(output.final$beta.chosen),"\n", sep = ""))}
  }
  return(output.final)
}

###################################################################################
## main function
###################################################################################
farm.select.adjusted <- function (Y, X,   loss=c("mcp", "scad", "lasso"),robust = FALSE, ...){
  n = length(Y)
  if(robust ==TRUE){
      Y.mean =mu_robust_F(matrix(Y,1,n), matrix(rep(1,n),n,1))
      Y=Y-rep(Y.mean,n)
      X.mean = mu_robust_F(matrix(X,p,n), matrix(rep(1,n),n,1))
      X = sweep(X, 1,X.mean,"-")
    }else{
      Y  = Y-mean(Y)
      X.mean = rowMeans(X)
      X = sweep(X, 1,X.mean,"-")
    }


  #adjust for factors

  output.adjust  = farm.adjust(Y=Y, X=t(X),robust= robust,...)
  #find residuals from factor adjustment
  X.res = output.adjust$X.res
  Y.res = output.adjust$Y.res
  nfactors = output.adjust$nfactors
  nx = NROW(X.res)
  p = NCOL(X.res)
  loss <- match.arg(loss)

  #perform model selection
  output.chosen = farm.select.temp (Y, X, X.res, Y.res, loss)

  #list all the output
  list(beta.chosen = output.chosen$beta_chosen, coef.chosen = output.chosen$coef_chosen, nfactors = nfactors, X.res = X.res, Y.res = Y.res)
}
#

# ################# adjusting factors ################
#' Adjusting a data matrix for underlying factors
#'
#' Given a matrix of covariates, this function estimates the underlying factors and computes data residuals after regressing out those factors.
#' @param Y a size n outcome vector, \emph{optional}.
#' @param X an n x p data matrix with each row being a sample. Must have same number of rows as the size of \code{Y}.
#' @param K.factors a \emph{optional} number of factors to be estimated. Otherwise estimated internally.
#' @param robust an \emph{optional} boolean, specifying whether or not to use robust estimators for mean and variance. Default is FALSE.
#' @return A list with the following items
#' \item{residual}{the data after being adjusted for underlying factors}
#' \item{loadings}{estimated factor loadings}
#' \item{factors}{estimated factors}
#' \item{nfactors}{the number of (estimated) factors}
#' \item{X.res}{the decorrelated covariate matrix}
#' \item{Y.res}{the response vector corresponding to the decorrelated covariates}
#' @details For details about the method, see Fan et al.(2017).
#' @details For formula of how the covariates and the repsonse vector are  adjusted for latent factors, see Section 3.2 in Fan et al.(2017).
# #' @details Using \code{robust.cov = TRUE} uses the Huber's loss to estimate the covariance matrix. For details of covariance estimation method see Fan et al.(2017).
#' @details Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details Number of latent factors, if not provided, is estimated by the eignevalue ratio test. See Ahn and Horenstein(2013).
#' @examples
#' set.seed(100)
#' P = 100 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 5 #model size
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' beta_1 = 3+3*runif(Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps = rnorm(N)
#' Y = X%*%beta+eps
#' output = farm.adjust(Y,X)
#'
#' @references Ahn, S. C., and A. R. Horenstein (2013): "Eigenvalue Ratio Test for the Number of Factors," Econometrica, 81 (3), 1203–1227.
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/abs/1612.08490}
#' @export
farm.adjust<- function(Y = NULL, X ,K.factors = NULL, robust = FALSE) {#, robust.cov = FALSE) {
  X = t(X)
  p = NROW(X)
  n = NCOL(X)
  if(min(n,p)<=4) stop('n and p must be at least 4')

  if(any(abs(rowMeans(X))>0.00000001)){
    if(robust ==TRUE){
      X.mean = mu_robust_F(matrix(X,p,n), matrix(rep(1,n),n,1))
      X = sweep(X, 1,X.mean,"-")
    }else{
      X.mean = rowMeans(X)
      X = sweep(X, 1,X.mean,"-")
    }
  }



  #estimate covariance matrix
 if(robust ==TRUE){
  covx =  Cov_Huber(matrix((X),p,n),  matrix(rep(1,n),n,1))
  }else{
  covx = cov(t(X))
  }
  eigs = Eigen_Decomp(covx)
  values = eigs[,p+1]
  vectors = eigs[,1:p]
  #estimate nfactors
  values = pmax(values,0)
  ratio=c()
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(n,p)/2))){
      ratio=append(ratio, values[i+1]/values[i])}
    ratio = ratio[is.finite(ratio)]
    K.factors = which.min(ratio)} else {K.factors}
  if(K.factors>min(n,p)/2) warning('Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies')
  if(K.factors>max(n,p)) stop('Number of factors cannot be larger than n or p')

  #using K.factors estimate the factors and loadings
  P_F = Find_factors( Sigma = matrix(covx,p,p), (X), n, p,  K.factors)
  X.res = Find_X_star (P_F, X)

  #output
  if(is.null(Y)){
  list(X.res = X.res, nfactors = K.factors)
  }else {
    if(NCOL(X)!=NROW(Y)) stop('number of rows in covariate matrix should be size of the outcome vector')
    if(abs(mean(Y))>0.00000001){
      if(robust ==TRUE){
        Y.mean =mu_robust_F(matrix(Y,1,n), matrix(rep(1,n),n,1))
        Y=Y-rep(Y.mean,n)
      }else{
        Y  = Y-mean(Y)
      }
    }
    Y.res = Find_Y_star ( P_F, matrix(Y,n,1))
    list(X.res = X.res, nfactors = K.factors, Y.res = Y.res)
  }
}

farm.select.temp<- function(Y, X, X.res, Y.res,loss)
{
  N = length(Y.res)
  if (loss == "scad"){
    CV_SCAD=cv.ncvreg(X.res, Y.res,penalty="SCAD",seed = 100, nfolds = ceiling(N/3))
    beta_SCAD=coef(CV_SCAD, s = "lambda.min", exact=TRUE)
    inds_SCAD=which(beta_SCAD!=0)
    inds_SCAD = inds_SCAD[ - which(inds_SCAD ==1)]
    coef_chosen = beta_SCAD[inds_SCAD]
    inds_SCAD = inds_SCAD-1
    list(beta_chosen= inds_SCAD, coef_chosen=coef_chosen )
  }
  else if (loss == "lasso"){
    CV_lasso=cv.ncvreg(X.res, Y.res,penalty="lasso", seed = 100,nfolds = ceiling(N/3))
    beta_lasso=coef(CV_lasso, s = "lambda.min", exact=TRUE)
    inds_lasso=which(beta_lasso!=0)
    inds_lasso = inds_lasso[ - which(inds_lasso ==1)]
    coef_chosen = beta_lasso[inds_lasso]
    inds_lasso = inds_lasso-1
    list(beta_chosen= inds_lasso, coef_chosen=coef_chosen )
  }else {

    CV_MCP=cv.ncvreg(X.res, Y.res,penalty="MCP",seed=100,nfolds = ceiling(N/3))
    beta_MCP=coef(CV_MCP, s = "lambda.min", exact=TRUE)
    inds_MCP=which(beta_MCP!=0)
    inds_MCP = inds_MCP[ - which(inds_MCP ==1)]
    coef_chosen = beta_MCP[inds_MCP]
    inds_MCP = inds_MCP-1
    list(beta_chosen= inds_MCP, coef_chosen=coef_chosen )

  }

}



#################### huber mean calculation ##############################################
#' Mean estimation with Huber's loss function
#'
#' This function estimates mean of multivariate data using the Huber's loss. The tuning parameter is chosen by cross validation.
#' @param X a n x p data matrix with each row being a sample.

#' @return A list with the following items
#' \item{muhat}{the mean vector}
#' @examples
#' set.seed(100)
#' p = 20
#' n = 10
#' X = matrix(rnorm( p*n, 0,1), nrow = n)
#' muhat = farm.mean(X)
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @export
farm.mean <- function (X){
  X = t(X)
  p  = NROW(X)
  n = NCOL(X)
  muhat = mu_robust_F(matrix(X,p,n), matrix(rep(1,n),n,1))
  return(muhat)
}



#################### huber covariance calculation ##############################################
#' Covariance estimation with Huber's loss function
#'
#' This function estimates covariance of multivariate data using the Huber's loss. The tuning parameter is chosen by cross validation.
#' @param X an n x p data matrix with each row being a sample.

#' @return A list with the following items
#' \item{covhat}{the covariance matrix}
#' @examples
#' set.seed(100)
#' p = 20
#' n = 10
#' X = matrix(rnorm( p*n, 0,1), nrow = n)
#' covhat = farm.cov(X)
#'
#' @references Huber, P.J. (1964). "Robust Estimation of a Location Parameter." The Annals of Mathematical Statistics, 35, 73–101.
#' @export
farm.cov <- function (X){
  X = t(X)
  p  = NROW(X)
  n = NCOL(X)
  covhat = Cov_Huber(matrix((X),p,n),  matrix(rep(1,n),n,1))
  return(covhat)
}
