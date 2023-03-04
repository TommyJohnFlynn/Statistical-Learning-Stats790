# clear 
if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

# set-up
library(splines)
library(Matrix)
set.seed(13)

# main function 
truncpolyspline <- function(x, df, natural=FALSE) {
  
  # truncated case
  if (natural == FALSE) {                            # altered from the notes 
    knots <- quantile(x, seq(0, 1, length = df - 1)) 
    trunc_fun <- function(k) {(x>=k)*(x-k)^3}
    S <- sapply(knots[1:(df-2)], trunc_fun)
    S <- as(S, "CsparseMatrix")
    S <- cbind(x, x^2, S)
  }
  
  # natural case 
  else{
    knots <- quantile(x, seq(0, 1, length = df - 1)) 
    K <- knots[df-1]
    K1 <- knots[df-2]
    d.K1 <- ((x>=K1)*(x-K1)^3 - (x>=K)*(x-K)^3)/(K-K1) # ESL Equation 5.5
    trunc_fun <- function(k) {
      ((x>=k)*(x-k)^3 - (x>=K)*(x-K)^3)/(K-k) - d.K1   # ESL Equation 5.4
    }
    S <- sapply(knots[1:(df-3)], trunc_fun)
    S <- as(S, "CsparseMatrix")
    S <- cbind(x, S)
  }
  
  return(S)
}

# sample data 
xvec <- seq(0, 1, length = 101) 
trunc <- truncpolyspline(sort(xvec), df=7, natural=FALSE)
nat <- truncpolyspline(sort(xvec), df=7, natural=TRUE)

# plot 
par(mfrow=c(1,2))
matplot(scale(trunc), type='l', main='Truncated Basis', 
        xlab='x', ylab='Scaled Transformation')
matplot(scale(nat), type='l',  main='Natural Basis',
        xlab='x', ylab='Scaled Transformation')
