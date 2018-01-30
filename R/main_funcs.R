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
#' @param verbose a boolean determining whether or not to print output to the console. Default is TRUE.
#' @param robust a boolean, specifying whether or not to use robust estimators for mean and variance. Default is FALSE.
#' @param lin.reg a boolean, specifying whether or not to assume that we have a linear regression model (TRUE) or a logit model (FALSE) structure. Default is TRUE.
#' @param K.factors a \emph{optional} number of factors to be estimated. Otherwise estimated internally. K>0.
#' @param max.iter Maximum number of iterations across the regularization path. Default is 10000.
#' @param nfolds The number of cross-validation folds. Default is ceiling(samplesize/3).
#' @return A list with the following items
#' \item{model.size}{the size of the model}
#' \item{beta.chosen}{the indices of the covariates chosen in the model}
#' \item{coef.chosen}{the coefficients of the chosen covariates}
#' \item{X.residual}{the residual covariate matrix after adjusting for factors}
#' \item{nfactors}{number of (estimated) factors}
#' @details Number of rows and columns of the covariate matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For details about the method, see Fan et al.(2017).
#' @details For formula of how the covariates are  adjusted for latent factors, see Section 3.2 in Fan et al.(2017).
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
#'
#' ##with default options
#' output = farm.select(Y,X)
#' output$beta.chosen #variables selected
#' output$coef.chosen #coefficients of selected variables
#'
#' ##changing the loss function and inputting factors
#' output = farm.select(Y,X, loss = "scad", K.factors = 4)
#'
#' ##use a logistic regression model
#' Prob = 1/(1+exp(-X%*%beta))
#' Y = rbinom(N, 1, Prob)
#' output = farm.select(Y,X, lin.reg=FALSE, loss = "lasso")
#'
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/abs/1612.08490}
#' @export
farm.select <- function (Y, X,  loss=c("mcp", "scad", "lasso"), verbose = TRUE, robust = FALSE,  lin.reg = TRUE,K.factors= NULL,max.iter=10000, nfolds=ceiling(length(Y)/3)){
  loss=match.arg(loss)
  #error checking
  X = t(X)
  if(NCOL(X)!=NROW(Y)) stop('number of rows in covariate matrix should be size of the outcome vector')
    output.final = farm.select.adjusted(Y, X,   loss=match.arg(loss), robust=robust,lin.reg=lin.reg,K.factors=K.factors, max.iter = max.iter, nfolds = nfolds)
    if(verbose){output.final.call = match.call()
    cat("\n Call:\n")
    print(output.final.call)
    cat("\n Factor Adjusted Robust Model Selection \n")
    cat(paste("loss function used: ",  match.arg(loss), "\n",sep = ""))
    cat(paste("\np = ", NROW(X),", n = ", NCOL(X),"\n", sep = ""))
    cat(paste("factors found: ", output.final$nfactors, "\n", sep = ""))
    cat("size of model selected:\n")
    cat(paste(" ", output.final$model.size,"\n", sep = ""))
  }
  return(output.final)
}

###################################################################################
## main function: not for end user
###################################################################################
farm.select.adjusted <- function (Y, X,   loss,robust ,lin.reg,K.factors, max.iter, nfolds){
  n = length(Y)
  p = NROW(X)
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

  output.adjust  = farm.adjust(Y=Y, X=t(X),robust= robust,lin.reg=lin.reg,K.factors)
  #find residuals from factor adjustment
  X.res = output.adjust$X.res
  Y.res = output.adjust$Y.res
  nfactors = output.adjust$nfactors
  nx = NROW(X.res)
  p = NCOL(X.res)
  X.residual = output.adjust$X.residual
  #perform model selection
  if(lin.reg){
   output.chosen = farm.select.temp ( X.res, Y.res, loss,max.iter, nfolds)
  }else{
    F_hat =  output.adjust$F_hat
    output.chosen = farm.select.temp ( X.res, Y.res, loss,max.iter, nfolds, factors = F_hat)
  }

  #list all the output
  list(model.size = output.chosen$model_size, beta.chosen = output.chosen$beta_chosen, coef.chosen = output.chosen$coef_chosen, nfactors = nfactors, X.residual =X.residual)
}
#

##adjusting for factors: internal function
farm.adjust<- function(Y , X , robust , lin.reg,K.factors ) {#, robust.cov = FALSE) {
  X = t(X)
  p = NROW(X)
  n = NCOL(X)
  if(min(n,p)<=4) stop('\n n and p must be at least 4 \n')

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
  if(K.factors>min(n,p)/2) warning('\n Warning: Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies \n')
  if(K.factors>max(n,p)) stop('\n Number of factors cannot be larger than n or p \n')

  #using K.factors estimate the factors and loadings
  Lambda_hat = Find_lambda_class(Sigma = matrix(covx,p,p), (X), n, p,  K.factors)
  F_hat = Find_factors_class(Lambda_hat, (X), n, p,  K.factors)
  X.res2 = Find_X_star_class(F_hat,Lambda_hat, X )

  if(lin.reg){
    P_F = Find_PF( F_hat,  n)
    X.res1 = Find_X_star( P_F, X)
  }


    if(NCOL(X)!=NROW(Y)) stop('\n number of rows in covariate matrix should be size of the outcome vector \n')
    if(abs(mean(Y))>0.00000001){
      if(robust ==TRUE){
        Y.mean =mu_robust_F(matrix(Y,1,n), matrix(rep(1,n),n,1))
        Y=Y-rep(Y.mean,n)
      }else{
        Y  = Y-mean(Y)
      }
    }
    if(lin.reg){
      Y.res1 = Find_Y_star( P_F, matrix(Y,n,1))
    }else{
      Y.res2 = matrix(Y,n,1)
    }
    if(lin.reg){
      list(X.res = X.res1,Y.res = Y.res1,  nfactors = K.factors, X.residual = X.res2)
    }else{
      list(X.res = X.res2,Y.res = Y.res2,  nfactors = K.factors, F_hat = F_hat,X.residual = X.res2)
    }
}
#model selection: not for end user
farm.select.temp<- function(X.res, Y.res,loss, max.iter, nfolds ,factors = NULL)
{

  N = length(Y.res)
  P = NCOL(X.res)
  if(is.null(factors)){
    if (loss == "scad"){
      CV_SCAD=cv.ncvreg(X.res, Y.res,penalty="SCAD",seed = 100, nfolds = nfolds, max.iter  = max.iter)
      beta_SCAD=coef(CV_SCAD, s = "lambda.min", exact=TRUE)
      inds_SCAD=which(beta_SCAD!=0)
      if( length(inds_SCAD)==1){
        list(beta_chosen= NULL, coef_chosen=NULL,model_size=0 )
        }else{
          inds_SCAD = inds_SCAD[ - which(inds_SCAD ==1)]
          coef_chosen = beta_SCAD[inds_SCAD]
          inds_SCAD = inds_SCAD-1
          list(beta_chosen= inds_SCAD, coef_chosen=coef_chosen,model_size=length(inds_SCAD))
        }
      }
    else if (loss == "lasso"){
      CV_lasso=cv.ncvreg(X.res, Y.res,penalty="lasso", seed = 100,nfolds = nfolds, max.iter  = max.iter)
      beta_lasso=coef(CV_lasso, s = "lambda.min", exact=TRUE)
      inds_lasso=which(beta_lasso!=0)
      if( length(inds_lasso)==1){
        list(beta_chosen= NULL, coef_chosen=NULL,model_size=0 )
        }else{
          inds_lasso = inds_lasso[ - which(inds_lasso ==1)]
          coef_chosen = beta_lasso[inds_lasso]
          inds_lasso = inds_lasso-1
          list(beta_chosen= inds_lasso, coef_chosen=coef_chosen,model_size=length(inds_lasso))
          }
      }else{
        CV_MCP=cv.ncvreg(X.res, Y.res,penalty="MCP",seed=100,nfolds = nfolds, max.iter  = max.iter)
        beta_MCP=coef(CV_MCP, s = "lambda.min", exact=TRUE)
        inds_MCP=which(beta_MCP!=0)
        if( length(inds_MCP)==1){
          list(beta_chosen= NULL, coef_chosen=NULL ,model_size=0)
          }else{
            inds_MCP = inds_MCP[ - which(inds_MCP ==1)]
            coef_chosen = beta_MCP[inds_MCP]
            inds_MCP = inds_MCP-1
            list(beta_chosen= inds_MCP, coef_chosen=coef_chosen ,model_size=length(inds_MCP))
          }
      }
    }else{
      Kfactors = NCOL(factors)
      if (loss == "scad"){
        CV_SCAD=cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="SCAD",seed = 100, nfolds = nfolds, family = "binomial", max.iter  = max.iter,penalty.factor = c(rep(0,Kfactors), rep(1, P)))
        beta_SCAD=coef(CV_SCAD, s = "lambda.min", exact=TRUE)
        beta_SCAD = beta_SCAD[-c(2:(Kfactors+1))]
        inds_SCAD=which(beta_SCAD!=0)
        if( length(inds_SCAD)==1){
         list(beta_chosen= NULL, coef_chosen=NULL,model_size=0 )
          }else{
            inds_SCAD = inds_SCAD[ - which(inds_SCAD ==1)]
            coef_chosen = beta_SCAD[inds_SCAD]
            inds_SCAD = inds_SCAD-1
            list(beta_chosen= inds_SCAD, coef_chosen=coef_chosen,model_size=length(inds_SCAD) )
          }
        }
      else if (loss == "lasso"){
        CV_lasso=cv.ncvreg(X = cbind(factors,X.res), y= Y.res,penalty="lasso", seed = 100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,penalty.factor = c(rep(0,Kfactors), rep(1, P)))
         beta_lasso=coef(CV_lasso, s = "lambda.min", exact=TRUE)
        beta_lasso = beta_lasso[-c(2:(Kfactors+1))]
        inds_lasso=which(beta_lasso!=0)
        if( length(inds_lasso)==1){
          list(beta_chosen= NULL, coef_chosen=NULL,model_size=0)
          }else{
            inds_lasso = inds_lasso[ - which(inds_lasso ==1)]
            coef_chosen = beta_lasso[inds_lasso]
            inds_lasso = inds_lasso-1
            list(beta_chosen= inds_lasso, coef_chosen=coef_chosen,model_size=length(inds_lasso) )
            }
        }else {
          CV_MCP=cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="MCP",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,penalty.factor = c(rep(0,Kfactors), rep(1, P)))
          beta_MCP=coef(CV_MCP, s = "lambda.min", exact=TRUE)
          beta_MCP = beta_MCP[-c(2:(Kfactors+1))]
          inds_MCP=which(beta_MCP!=0)
          if( length(inds_MCP)==1){
            list(beta_chosen= NULL, coef_chosen=NULL, model_size =0 )
            }else{
              inds_MCP = inds_MCP[ - which(inds_MCP ==1)]
              coef_chosen = beta_MCP[inds_MCP]
              inds_MCP = inds_MCP-1
              list(beta_chosen= inds_MCP, coef_chosen=coef_chosen, model_size =length(inds_MCP) )
            }
        }
    }
  }



# ################# adjusting factors ################
#' Adjusting a data matrix for underlying factors
#'
#' Given a matrix of covariates, this function estimates the underlying factors and computes data residuals after regressing out those factors.
#' @param X an n x p data matrix with each row being a sample.
#' @param K.factors a \emph{optional} number of factors to be estimated. Otherwise estimated internally. K>0.
#' @param robust an \emph{optional} boolean, specifying whether or not to use robust estimators for mean and variance. Default is FALSE.
#' @return A list with the following items
#' \item{residual}{the data after being adjusted for underlying factors}
#' \item{loadings}{estimated factor loadings}
#' \item{factors}{estimated factors}
#' \item{nfactors}{the number of (estimated) factors}
#' @details For details about the method, see Fan et al.(2017).
#' @details Using \code{robust.cov = TRUE} uses the Huber's loss to estimate the covariance matrix. For details of covariance estimation method see Fan et al.(2017).
#' @details Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details Number of latent factors, if not provided, is estimated by the eignevalue ratio test. See Ahn and Horenstein(2013).
#' @examples
#' set.seed(100)
#' P = 100 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' output = farm.res(X) #default options
#' output$nfactors
#' output = farm.res(X, K.factors = 10) #inputting factors
#' names(output) #list of output
#' @references Ahn, S. C., and A. R. Horenstein (2013): "Eigenvalue Ratio Test for the Number of Factors," Econometrica, 81 (3), 1203–1227.
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/abs/1612.08490}
#' @export
farm.res<- function(X ,K.factors = NULL, robust = FALSE) {#, robust.cov = FALSE) {
  X = t(X)
  p = NROW(X)
  n = NCOL(X)
  if(min(n,p)<=4) stop('\n n and p must be at least 4 \n')

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
  if(K.factors>min(n,p)/2) warning('\n Warning: Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies \n')
  if(K.factors>max(n,p)) stop('\n Number of factors cannot be larger than n or p \n')

  #using K.factors estimate the factors and loadings
  Lambda_hat = Find_lambda_class(Sigma = matrix(covx,p,p), (X), n, p,  K.factors)
  F_hat = Find_factors_class(Lambda_hat, (X), n, p,  K.factors)
  X.res = Find_X_star_class(F_hat,Lambda_hat, X )

  list(X.res = X.res, nfactors = K.factors, factors = F_hat, loadings  = Lambda_hat)

}
