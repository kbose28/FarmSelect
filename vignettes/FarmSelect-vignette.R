## ----gh-installation, eval = FALSE---------------------------------------
#  install.packages("devtools")
#  devtools::install_github("kbose28/FarmSelect")
#  library(FarmSelect)

## ------------------------------------------------------------------------
library(FarmSelect)
set.seed(100)
P = 200 #dimension
N = 50 #samples
K = 3 #nfactors
Q = 5 #model size
Lambda = matrix(rnorm(P*K, 0,1), P,K)
F = matrix(rnorm(N*K, 0,1), N,K)
U = matrix(rnorm(P*N, 0,1), P,N)
X = Lambda%*%t(F)+U
X = t(X)
beta_1 = 3+3*runif(Q)
beta = c(beta_1, rep(0,P-Q))
eps = rnorm(N)
Y = X%*%beta+eps 
output = farm.select(X,Y)

## ------------------------------------------------------------------------
names(output)
output$beta.chosen
output$coef.chosen

## ------------------------------------------------------------------------
output = farm.select(X,Y, loss = "mcp", K.factors = 10)

## ------------------------------------------------------------------------
set.seed(100)
Prob = 1/(1+exp(-X%*%beta))
Y = rbinom(N, 1, Prob)
output = farm.select(X,Y, lin.reg=FALSE, loss = "lasso")

## ------------------------------------------------------------------------
output = farm.res(X)
names(output)

## ------------------------------------------------------------------------
output = farm.res(X, K.factors  = 30)

