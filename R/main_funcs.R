#' @useDynLib FarmSelect
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics mtext plot points axis par barplot
#' @importFrom fBasics rowSds
#' @import methods
#' @import utils
#' @import grDevices
#' @import stats
#' @import ncvreg

NULL
###################################################################################
## Main function that performs model selection
###################################################################################
#' Factor-adjusted robust model selection
#'
#' Given a covariate matrix and output vector, this function first adjusts the covariates for underlying factors and then performs model selection.
#' @param X an n x p covariate matrix with each row being a sample. Must have same number of rows as the size of \code{Y}.
#' @param Y a size n outcome vector.
#' @param loss a character string specifying the loss function to be minimized. Must be one of "scad" (default) "mcp" or "lasso". You can just specify the initial letter.
#' @param robust a boolean, specifying whether or not to use robust estimators for mean and variance. Default is TRUE.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is FALSE. Only used if \code{robust} is TRUE.
#' @param tau \code{>0}, multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{robust} is TRUE and \code{cv} is FALSE. See details.
#' @param lin.reg a boolean, specifying whether or not to assume that we have a linear regression model (TRUE) or a logit model (FALSE) structure. Default is TRUE.
#' @param K.factors number of factors to be estimated. Otherwise estimated internally. K>0.
#' @param max.iter maximum number of iterations across the regularization path. Default is 10000.
#' @param nfolds the number of cross-validation folds. Default is ceiling(samplesize/3).
#' @param eps Convergence threshhold for model fitting using \code{\link{ncvreg}}. The algorithm iterates until the RMSD for the change in linear predictors for any coefficient is less than eps. Default is \code{1e-4}.
#' @param verbose a boolean specifying whether to print runtime updates to the console. Default is TRUE.
#' @return A list with the following items
#' \item{model.size}{the size of the model}
#' \item{beta.chosen}{the indices of the covariates chosen in the model}
#' \item{coef.chosen}{the coefficients of the chosen covariates}
#' \item{X.residual}{the residual covariate matrix after adjusting for factors}
#' \item{nfactors}{number of (estimated) factors}
#' \item{n}{number of observations}
#' \item{p}{number of dimensions}
#' \item{robust}{whether robust parameters were used}
#' \item{loss}{loss function used}
#' #' @details Number of rows and columns of the covariate matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For formula of how the covariates are  adjusted for latent factors, see Section 3.2 in Fan et al.(2017).
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate } is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @details \code{\link{ncvreg}} is used to fit the model after decorrelation. This package may output its own warnings about failures to converge and model saturation.
#' @examples
#' ##linear regression
#' set.seed(100)
#' P = 200 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 3 #model size
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' beta_1 = rep(5,Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps = rt(N, 2.5)
#' Y = X%*%beta+eps
#'
#' ##with default options
#' output = farm.select(X,Y) #robust, no cross-validation
#' output$beta.chosen #variables selected
#' output$coef.chosen #coefficients of selected variables
#'
#' #examples of other robustification options
#' output = farm.select(X,Y,robust = FALSE) #non-robust
#' output = farm.select(X,Y, tau = 3) #robust, no cross-validation, specified tau
#' #output = farm.select(X,Y, cv= TRUE) #robust, cross-validation: LONG RUNNING!
#'
#' ##changing the loss function and inputting factors
#' output = farm.select(X, Y,loss = "mcp", K.factors = 4)
#'
#' ##use a logistic regression model, a larger sample size is desired.
#' \dontrun{
#' set.seed(100)
#' P = 400 #dimension
#' N = 300 #samples
#' K = 3 #nfactors
#' Q = 3 #model size
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' beta_1 = rep(5, Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps = rnorm(N)
#' Prob = 1/(1+exp(-X%*%beta))
#' Y = rbinom(N, 1, Prob)
#'
#' output = farm.select(X,Y, lin.reg=FALSE, eps=1e-3)
#' output$beta.chosen
#' output$coef.chosen
#' }
#' @seealso \code{\link{print.farm.select}} \code{\link{farm.res}}
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/abs/1612.08490}
#' @export
farm.select <- function ( X, Y, loss=c( "scad","mcp", "lasso"),  robust = TRUE, cv = FALSE, tau = 2 ,  lin.reg = TRUE,K.factors= NULL,max.iter=10000, nfolds=ceiling(length(Y)/3), eps =1e-4,verbose=TRUE){
  loss=match.arg(loss)
  #error checking
  X = t(X)
  if(tau<=0) stop('tau should be a positive number')
  if(NCOL(X)!=NROW(Y)) stop('number of rows in covariate matrix should be size of the outcome vector')
  output.final = farm.select.adjusted( X, Y,  loss=match.arg(loss), robust=robust,cv=cv, tau=tau, lin.reg=lin.reg,K.factors=K.factors, max.iter = max.iter, nfolds = nfolds,eps=eps, verbose=verbose)
  output = (output.final)
  attr(output, "class") <- "farm.select"
  output
}
###################################################################################
## Print method
###################################################################################
#' Summarize and print the results of the model selection
#'
#' Print method for \code{farm.select} objects
#' @param x A \code{farm.select} object.
#' @param \dots Further arguments passed to or from other methods.
#' @return A list with the following items:
#' \item{model.size}{the size of the model}
#' \item{beta.chosen}{the indices of the covariates chosen in the model}
#' \item{coef.chosen}{the coefficients of the chosen covariates}
#' \item{X.residual}{the residual covariate matrix after adjusting for factors}
#' \item{nfactors}{number of (estimated) factors}
#' \item{n}{number of observations}
#' \item{p}{number of dimensions}
#' \item{robust}{whether robust parameters were used}
#' \item{loss}{loss function used}
#' @seealso \code{\link{farm.select}}
#'
#' @examples
#' set.seed(100)
#' P = 200 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 3 #model size
#' Lambda = matrix(rnorm(P*K, 0,1), P,K)
#' F = matrix(rnorm(N*K, 0,1), N,K)
#' U = matrix(rnorm(P*N, 0,1), P,N)
#' X = Lambda%*%t(F)+U
#' X = t(X)
#' beta_1 = rep(5,Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps = rt(N, 2.5)
#' Y = X%*%beta+eps
#'
#' ##with default options
#' output = farm.select(X,Y)
#' output
#' @export
print.farm.select<-function(x,...){
  cat(paste("\n Factor Adjusted",if(x$robust) "Robust", "Model Selection \n"))
  cat(paste("loss function used: ",  x$loss, "\n",sep = ""))
  cat(paste("\np = ", x$p,", n = ", x$n,"\n", sep = ""))
  cat(paste("factors found: ", x$nfactors, "\n", sep = ""))
  cat("size of model selected:\n")
  cat(paste(" ", x$model.size,"\n", sep = ""))

 }


###################################################################################
## main function: not for end user
###################################################################################
farm.select.adjusted <- function (X,  Y,  loss,robust ,cv=cv,tau, lin.reg,K.factors, max.iter, nfolds,eps,verbose){
  n = length(Y)
  p = NROW(X)

  if(robust ==TRUE){
    if(cv==TRUE){
     if(lin.reg==TRUE){Y.mean =mu_robust_F(matrix(Y,1,n), matrix(rep(1,n),n,1))
        Y=Y-rep(Y.mean,n)}
        X.mean = mu_robust_F(matrix(X,p,n), matrix(rep(1,n),n,1))
        X = sweep(X, 1,X.mean,"-")
        }else{
         if(lin.reg==TRUE){
            CT = tau*sd(Y)*sqrt(n/log(n))
            Y.mean = mu_robust_F_noCV(matrix(Y,1,n), matrix(rep(1,n),n,1), matrix(CT,1,1))
            Y=Y-rep(Y.mean,n)}
            CT = tau*rowSds(X)*sqrt(n/log(p*n))
            X.mean = mu_robust_F_noCV(matrix(X,p,n), matrix(rep(1,n),n,1), matrix(CT, p,1))
            X = sweep(X, 1,X.mean,"-")}
  }else{
    if(lin.reg==TRUE){Y  = Y-mean(Y)}
    X.mean = rowMeans(X)
    X = sweep(X, 1,X.mean,"-")
  }


  #adjust for factors

  output.adjust  = farm.adjust(Y=Y, X=t(X),robust= robust,cv=cv,tau = tau, lin.reg=lin.reg,K.factors, verbose)
  #find residuals from factor adjustment
  X.res = output.adjust$X.res
  Y.res = output.adjust$Y.res
  nfactors = output.adjust$nfactors
  nx = NROW(X.res)
  p = NCOL(X.res)
  X.residual = output.adjust$X.residual
  #perform model selection
  if(lin.reg){
    output.chosen = farm.select.temp ( X.res,Y.res,  loss,max.iter, nfolds, eps)
  }else{
    F_hat =  output.adjust$F_hat
    output.chosen = farm.select.temp (  X.res,Y.res, loss,max.iter, nfolds, eps,factors = F_hat)
  }

  #list all the output
  list(model.size = output.chosen$model_size, beta.chosen = unname(output.chosen$beta_chosen), coef.chosen = unname(output.chosen$coef_chosen), nfactors = nfactors, X.residual =unname(X.residual), p = p, n=n, robust = robust, loss =loss)
}
#

##adjusting for factors: internal function
farm.adjust<- function( X ,Y , robust ,cv=cv,tau  = tau,  lin.reg,K.factors, verbose ) {
  X = t(X)
  p = NROW(X)
  n = NCOL(X)
  if(min(n,p)<=4) stop('\n n and p must be at least 4 \n')


  #estimate covariance matrix
  if(robust ==TRUE){
    if(cv==TRUE){if(verbose){cat("calculating covariance matrix...\n")}
    covx =  Cov_Huber(matrix((X),p,n),  matrix(rep(1,n),n,1))
    eigs = Eigen_Decomp(covx)
    values = eigs[,p+1]
    vectors = eigs[,1:p]
    }else{ if (verbose){cat("calculating tuning parameters...\n")}
        CT = Cov_Huber_tune(X, tau)
       if (verbose){cat("calculating covariance matrix...\n")}
        covx =  Cov_Huber_noCV(matrix((X),p,n),  matrix(rep(1,n),n,1), matrix(CT,p,p))
      eigs = Eigen_Decomp(covx)
      values = eigs[,p+1]
      vectors = eigs[,1:p]
    }
  }else{
    covx = cov(t(X))
    pca_fit=stats::prcomp((X), center = TRUE, scale = TRUE)
    values = (pca_fit$sdev)^2
    vectors = pca_fit$rotation
  }
if(verbose) { cat("fitting model...\n")}

  #estimate nfactors
  values = pmax(values,0)
  ratio=c()
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(n,p)/2))){
      ratio=append(ratio, values[i+1]/values[i])}
    ratio = ratio[is.finite(ratio)]
    K.factors = which.min(ratio)} else {K.factors}
  if(K.factors>min(n,p)/2) warning('\n Warning: Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies \n')
  if(K.factors>max(n,p)) stop('\n Number of factors supplied cannot be larger than n or p \n')

  #using K.factors estimate the factors and loadings
  Lambda_hat = Find_lambda_class(Sigma = matrix(covx,p,p), (X), n, p,  K.factors)
  F_hat = Find_factors_class(Lambda_hat, (X), n, p,  K.factors)
  X.res2 = Find_X_star_class(F_hat,Lambda_hat, X )

  if(lin.reg){
    P_F = Find_PF( F_hat,  n)
    X.res1 = Find_X_star( P_F, X)
  }


  if(NCOL(X)!=NROW(Y)) stop('\n number of rows in covariate matrix should be size of the outcome vector \n')

  if(lin.reg){
    Y.res1 = Find_Y_star( P_F, matrix(Y,n,1))
  }else{
    Y.res2 = matrix(Y,n,1)
  }
  if(lin.reg){
    list( X.res = X.res1,Y.res = Y.res1, nfactors = K.factors, X.residual = X.res2)
  }else{
    list(X.res = X.res2, Y.res = Y.res2, nfactors = K.factors, F_hat = F_hat,X.residual = X.res2)
  }
}
#model selection: not for end user
farm.select.temp<- function( X.res,Y.res,loss, max.iter, nfolds ,eps,factors = NULL)
{

  N = length(Y.res)
  P = NCOL(X.res)
  if(is.null(factors)){
    if (loss == "mcp"){
      CV_SCAD=cv.ncvreg(X.res, Y.res,penalty="MCP",seed = 100, nfolds = nfolds, max.iter  = max.iter, eps=eps)
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
      CV_lasso=cv.ncvreg(X.res, Y.res,penalty="lasso", seed = 100,nfolds = nfolds, max.iter  = max.iter,eps=eps)
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
      CV_MCP=cv.ncvreg(X.res, Y.res,penalty="SCAD",seed=100,nfolds = nfolds, max.iter  = max.iter,eps=eps)
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
    if (loss == "mcp"){
      CV_MCP=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="MCP",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps,penalty.factor = c(rep(0,Kfactors), rep(1, P))), silent=TRUE)
      if("try-error" %in% class(t)){CV_MCP=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="MCP",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps))}
      beta_MCP=coef(CV_MCP, s = "lambda.min", exact=TRUE)
      beta_MCP = beta_MCP[-c(2:(Kfactors+1))]
      inds_MCP=which(beta_MCP!=0)
      if( length(inds_MCP)==1){
        list(beta_chosen= NULL, coef_chosen=NULL,model_size=0 )
      }else{
        inds_MCP = inds_MCP[ - which(inds_MCP ==1)]
        coef_chosen = beta_MCP[inds_MCP]
        inds_MCP = inds_MCP-1
        list(beta_chosen= inds_MCP, coef_chosen=coef_chosen,model_size=length(inds_MCP) )
      }
    }
    else if (loss == "lasso"){
      CV_lasso=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="lasso",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps,penalty.factor = c(rep(0,Kfactors), rep(1, P))), silent=TRUE)
      if (inherits(CV_lasso, "try-error")){CV_lasso=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="lasso",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps))}
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
      CV_SCAD=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="SCAD",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps,penalty.factor = c(rep(0,Kfactors), rep(1, P))), silent=TRUE)
      if("try-error" %in% class(t)){CV_SCAD=try(cv.ncvreg(X = cbind(factors,X.res), y = Y.res,penalty="SCAD",seed=100,nfolds = nfolds, family = "binomial", max.iter  = max.iter,eps=eps))}
      beta_SCAD=coef(CV_SCAD, s = "lambda.min", exact=TRUE)
      beta_SCAD = beta_SCAD[-c(2:(Kfactors+1))]
      inds_SCAD=which(beta_SCAD!=0)
      if( length(inds_SCAD)==1){
        list(beta_chosen= NULL, coef_chosen=NULL, model_size =0 )
      }else{
        inds_SCAD = inds_SCAD[ - which(inds_SCAD ==1)]
        coef_chosen = beta_SCAD[inds_SCAD]
        inds_SCAD = inds_SCAD-1
        list(beta_chosen= inds_SCAD, coef_chosen=coef_chosen, model_size =length(inds_SCAD) )
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
#' @param robust a boolean, specifying whether or not to use robust estimators for mean and variance. Default is TRUE.
#' @param cv a boolean, specifying whether or  not to run cross-validation for the tuning parameter. Default is FALSE. Only used if \code{robust} is TRUE.
#' @param tau \code{>0} multiplier for the tuning parameter for Huber loss function. Default is 2. Only used if \code{robust} is TRUE and \code{cv} is FALSE. See details.
#' @param verbose a boolean specifying whether to print runtime updates to the console. Default is TRUE.
#' @return A list with the following items
#' \item{residual}{the data after being adjusted for underlying factors}
#' \item{loadings}{estimated factor loadings}
#' \item{factors}{estimated factors}
#' \item{nfactors}{the number of (estimated) factors}
#' @details For details about the method, see Fan et al.(2017).
#' @details Using \code{robust = TRUE} uses the Huber's loss to estimate parameters robustly. For details of covariance estimation method see Fan et al.(2017).
#' @details Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details Number of latent factors, if not provided, is estimated by the eignevalue ratio test. See Ahn and Horenstein(2013). The maximum number is taken to be min(n,p)/2. User can supply a larger number is desired.
#' @details The tuning parameter \code{= tau *  sigma * optimal rate } where \code{optimal rate} is the optimal rate for the tuning parameter. For details, see Fan et al.(2017). \code{sigma} is the standard deviation of the data.
#' @seealso \code{\link{farm.select}}
#' @examples
#' set.seed(100)
#' P = 200 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 3 #model size
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
farm.res<- function(X ,K.factors = NULL, robust = TRUE, cv=FALSE, tau = 2, verbose = TRUE) {
  X = t(X)
  p = NROW(X)
  n = NCOL(X)
  if(min(n,p)<=4) stop('\n n and p must be at least 4 \n')
if(tau<=0) stop('tau should be a positive number')
  if(robust ==TRUE){
    if(cv==TRUE){
    X.mean = mu_robust_F_noCV(matrix(X,p,n), matrix(rep(1,n),n,1))
    X = sweep(X, 1,X.mean,"-")
    }else{
       CT = tau*rowSds(X)*max(sqrt(n/log(p*n)))
      X.mean = mu_robust_F_noCV(matrix(X,p,n), matrix(rep(1,n),n,1), matrix(CT, p,1))
      X = sweep(X, 1,X.mean,"-")
    }
  }else{X.mean = rowMeans(X)
  X = sweep(X, 1,X.mean,"-")
  }

  #estimate covariance matrix
  if(robust ==TRUE){

 if(cv==TRUE){
   if(verbose){cat("calculating covariance matrix...\n")}
   covx =  Cov_Huber_noCV(matrix((X),p,n),  matrix(rep(1,n),n,1) )
   eigs = Eigen_Decomp(covx)
   values = eigs[,p+1]
   vectors = eigs[,1:p]
   } else{if(verbose){cat("calculating tuning parameters...\n")}
    CT = Cov_Huber_tune(X, tau)
   if(verbose){ cat("calculating covariance matrix...\n")}
    covx =  Cov_Huber_noCV(matrix((X),p,n),  matrix(rep(1,n),n,1), matrix(CT,p,p))
    eigs = Eigen_Decomp(covx)
    values = eigs[,p+1]
    vectors = eigs[,1:p]}
  }else{
    covx = cov(t(X))
    pca_fit=stats::prcomp((X), center = TRUE, scale = TRUE)
    values = (pca_fit$sdev)^2
    vectors = pca_fit$rotation
  }


  #estimate nfactors
  values = pmax(values,0)
  ratio=c()
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(n,p)/2))){
      ratio=append(ratio, values[i+1]/values[i])}
    ratio = ratio[is.finite(ratio)]
    K.factors = which.min(ratio)} else {K.factors}
  if(K.factors>min(n,p)/2) warning('\n Warning: Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies \n')
  if(K.factors>max(n,p)) stop('\n Number of factors supplied cannot be larger than n or p \n')

  #using K.factors estimate the factors and loadings
  Lambda_hat = Find_lambda_class(Sigma = matrix(covx,p,p), (X), n, p,  K.factors)
  F_hat = Find_factors_class(Lambda_hat, (X), n, p,  K.factors)
  X.res = Find_X_star_class(F_hat,Lambda_hat, X )

  list(X.res = X.res, nfactors = K.factors, factors = F_hat, loadings  = Lambda_hat)

}
