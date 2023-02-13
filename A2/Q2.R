if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

# set-up
library(glmnet)
library(microbenchmark)
set.seed(13)

# data augmentation  
augment <- function(X, Y, lambda){
  p <- ncol(X)
  X.aug <- rbind(X, sqrt(lambda)*diag(p))
  Y.aug <- c(Y, rep(0, p))
  return(list(X=X.aug, Y=Y.aug))
}

# augmented fit 
aug.fit <- function(X, Y, lambda){
  D = augment(scale(X, center=TRUE, scale=FALSE), Y, lambda)
  return(lm.fit(D$X, D$Y)$coefficients)
}

# load data 
prostate <- read.delim('prostate.txt')
X <- as.matrix(prostate[,2:9])
Y <-  prostate$lpsa

# augmented vs native fits 
lambda1 <- 0.1
lambda2 <- 0.1/(4*nrow(X)) # looks like they used different formula 
aug <- aug.fit(X,Y, lambda1)
native <- glmnet(X, Y, family='gaussian', alpha=0, lambda=lambda2)

# results 
as.vector(aug)         # intercept = 0.181560845
as.vector(native$beta) # intercept = 0.178745157 


# test run-times on large matrix 
generate.data <- function(n, p){
  Y <- rnorm(n)
  X <- matrix(rnorm(p*n), ncol = p)
  list(X = X, Y = Y)
}

D5 <- generate.data(10^5, 10)
m <- microbenchmark(aug.fit(D5$X,D5$Y, lambda1), 
               glmnet(D5$X, D5$Y, family='gaussian', alpha=0, lambda=lambda2),
               times= 100,
               unit='ms')
summary(m)$mean
