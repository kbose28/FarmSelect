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
#' @param Y a n x 1 outcome vector.
#' @param X a p x n covariate matrix with each row being a sample. Must have same number of columns as the size of \code{Y}.
#' @param loss an \emph{optional} character string specifying the loss function to be minimized. Must be one of "lasso" (default), "scad" or "mcp". You can just specify the initial letter.
#' @param verbose an \emph{optional} logical determining whether or not to print output to the console. Default is TRUE.
#' @param \dots Arguments passed to the \code{\link{farm.adjust}} function.
#' @return A list with the following items
#' \item{beta.chosen}{the indices of the covariates chosen in the model}
#' \item{coef.chosen}{the coefficients of the chosen coavriates}
#' \item{nfactors}{number of (estimated) factors}
#' @details
# #' @details
#' Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details For details about the method, see Fan et al.(2017).
#'
#' @examples
#' set.seed(100)
#' P = 100 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 5 #model size
#' Lambda=matrix(rnorm(P*K, 0,1), P,K)
#' F=matrix(rnorm(N*K, 0,1), N,K)
#' UU=matrix(rnorm(P*N, 0,1), P,N)
#' X=Lambda%*%t(F)+UU
#' beta_1=3+3*runif(Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps=rnorm(N)
#' Y=t(X)%*%beta+eps
#' output = farm.select(Y,X)
#' @references Fan J., Ke Y., Wang K., "Decorrelation of Covariates for High Dimensional Sparse Regression." \url{https://arxiv.org/pdf/1612.08490.pdf}
#' @export
farm.select <- function (Y, X,  loss=c("lasso", "scad", "mcp"), verbose = TRUE, ...){
  #error checking
  if(NCOL(X)!=NROW(Y)) stop('number of rows in covariate matrix should be size of the outcome vector')
    output.final = farm.select.adjusted(Y, X,   loss=c("lasso", "scad", "mcp"), ...)
    if(verbose){output.final.call = match.call()
    cat("Call:\n")
    print(output.final.call)
    cat("\n Factor Adjusted Robust Model Selection \n")
    cat(paste("loss function used: ",  match.arg(loss), "\n",sep = ""))
    cat(paste("\np = ", NROW(X),", n = ", NCOL(X),"\n", sep = ""))
    cat("size of model selected:\n")
    if(is.character(output.final$beta.chosen)){ cat(" no variable selected\n")} else{ cat(paste(" ", NROW(output$beta.chosen),"\n", sep = ""))}
  }
  return(output)
}

###################################################################################
## main function
###################################################################################
farm.select.adjusted <- function (Y, X,   loss=c("lasso", "scad", "mcp"), ...){
  X = sweep(X, 1, apply(X, 1, mean))
  #adjust for factors
  output.adjust  = farm.adjust(X, Y,...)

  #find residuals from factor adjustment
  U = output.adjust$X.res
  Y_star = output.adjust$Y.res
  nfactors = output.adjust$nfactors
  nx = NROW(U)
  p = NCOL(U)
  loss <- match.arg(loss)

  #perform model selection
  output.chosen = farm.select.temp (X, Y, U, Y_star, loss)

  #list all the output
  list(beta.chosen = output.chosen$beta_chosen, coef.chosen = output.chosen$coef_chosen, nfactors = nfactors)
}
#

# ################# adjusting factors ################
#' Adjusting a data matrix for underlying factors
#'
#' Given a matrix of covariates, this function estimates the underlying factors and computes data residuals after regressing out those factors.
#' @param Y a n x 1 outcome vector.
#' @param X a n x p data matrix with each row being a sample.
#' @param K.factors a \emph{optional} number of factors to be estimated. Otherwise estimated internally.
# #' @param robust.cov an \emph{optional} boolean, specifying whether or not to use a robust covariance matrix for factor adjustment. Default is FALSE.
#' @return A list with the following items
#' \item{residual}{the data after being adjusted for underlying factors}
#' \item{loadings}{estimated factor loadings}
#' \item{factors}{estimated factors}
#' \item{nfactors}{the number of (estimated) factors}

# #' @details Using \code{robust.cov = TRUE} uses the Huber's loss to estimate the covariance matrix. For details of covariance estimation method see Fan et al.(2017).
#' @details Number of rows and columns of the data matrix must be at least 4 in order to be able to calculate latent factors.
#' @details Number of latent factors, if not provided, is estimated by the eignevalue ratio test. See Ahn and Horenstein(2013).
#' @examples
#' set.seed(100)
#' P = 100 #dimension
#' N = 50 #samples
#' K = 3 #nfactors
#' Q = 5 #model size
#' Lambda=matrix(rnorm(P*K, 0,1), P,K)
#' F=matrix(rnorm(N*K, 0,1), N,K)
#' UU=matrix(rnorm(P*N, 0,1), P,N)
#' X=Lambda%*%t(F)+UU
#' beta_1=3+3*runif(Q)
#' beta = c(beta_1, rep(0,P-Q))
#' eps=0.1*rnorm(N)
#' Y=t(X)%*%beta+eps
#' output = farm.adjust(X,Y)
#'
#' @references Ahn, S. C., and A. R. Horenstein (2013): "Eigenvalue Ratio Test for the Number of Factors," Econometrica, 81 (3), 1203–1227.
#' @export
farm.adjust<- function(Y, X ,K.factors = NULL) {#, robust.cov = FALSE) {

  nx = NCOL(X)
  p = NROW(X)
  if(min(nx,p)<=4) stop('n and p must be at least 4')

  #estimate covariance matrix
  covx = cov(t(X))
  eigs = Eigen_Decomp( covx)
  values = eigs[,p+1]
  vectors = eigs[,1:p]

  #estimate nfactors
  values = pmax(values,0)
  ratio=c()
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(nx,p)/2))){
      ratio=append(ratio, values[i+1]/values[i])}
    ratio = ratio[is.finite(ratio)]
    K.factors = which.min(ratio)} else {K.factors}
  if(K.factors>min(nx,p)/2) warning('Number of factors supplied is > min(n,p)/2. May cause numerical inconsistencies')
  if(K.factors>max(nx,p)) stop('Number of factors cannot be larger than n or p')

  #using K.factors estimate the factors and loadings
  P_F = Find_factors( (X), nx, p,  K.factors)
  Y_star = Find_Y_star ( P_F, Y)
  U = Find_loading (P_F, X)

  #output
  list(X.res = U, nfactors = K.factors, Y.res = Y_star)
}

farm.select.temp<- function(X, Y, U, Y_star,loss)
{

  if (loss == "scad"){
    CV_SCAD=cv.ncvreg(U, Y_star,penalty="SCAD")
    beta_SCAD=coef(CV_SCAD, s = "lambda.min", exact=TRUE)
    inds_SCAD=which(beta_SCAD!=0)
    inds_SCAD = inds_SCAD[ - which(inds_SCAD ==1)]
    inds_SCAD=inds_SCAD-1
    list(beta_chosen= inds_SCAD, coef_chosen=beta_SCAD[beta_SCAD!=0] )
  }
  else if (loss == "mcp"){
    CV_MCP=cv.ncvreg(U, Y_star,penalty="MCP")
    beta_MCP=coef(CV_MCP, s = "lambda.min", exact=TRUE)
    inds_MCP=which(beta_MCP!=0)
    inds_MCP = inds_MCP[ - which(inds_MCP ==1)]
    inds_MCP=inds_MCP-1
    list(beta_chosen= inds_MCP, coef_chosen=beta_MCP[beta_MCP!=0] )
  }else {
    CV_1=cv.glmnet(U, Y_star)
    beta_1=as.vector(coef(CV_1, s = "lambda.min", exact=TRUE))
    inds_1=which(beta_1!=0)
    inds_1 = inds_1[ - which(inds_1 ==1)]
    inds_1=inds_1-1
    list(beta_chosen= inds_1, coef_chosen=beta_1[beta_1!=0] )
  }


}
