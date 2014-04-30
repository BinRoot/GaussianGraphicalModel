library(R.matlab)
library(matrixStats)
library(rbenchmark)

source('problem_parameters.R')
source('data_parameters.R')
source('ADMM_CNJGL/ADM_CNJGL.R')
source('ADMM_CNJGL/expand_ni.R')
source('ADMM_CNJGL/soft_scal.R')
source('ADMM_CNJGL/ADM_parameters.R')
source('ADMM_CNJGL/nodiag_construction.R')
source('ADMM_CNJGL/soft_thresh.R')
source('ADMM_CNJGL/diag_construction.R')
source('ADMM_CNJGL/soft_scal_nodiag.R')

p <- 50

data <- readMat('./data/X.mat')
X <- data$X

data <- readMat('./data/Y.mat')
Y <- data$Y
Y <- Y[,2:2]

# all -1 value entires
X1 <- X[which(Y %in% 1),]


# all 1 value entires
X2 <- X[which(Y %in% 1),]


x1rows = dim(X1)[1]
x1cols = dim(X1)[2]

x2rows = dim(X2)[1]
x2cols = dim(X2)[2]

indices <- order(colSds(X))
X1 <- X1[1:x1rows, indices[(length(indices)-9):length(indices)]]
X2 <- X2[1:x2rows, indices[(length(indices)-9):length(indices)]]

S1 <- cov(X1)

p <- dim(S1)[1]


S2 <- cov(X2)

print("input is")
print(S1)

t <- benchmark({
  res <- ADM_CNJGL(S1, S2, lambda_1, lambda_2, n1, n2)
}, replications=0)

print(t)

