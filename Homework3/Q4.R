# clear 
if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

# set-up
library(mgcv)
library(microbenchmark)
set.seed(13)

# generates data from a spherical cap on the unit sphere - from top to pi/3
surface.data <- function(n) {
  
  # 2D parameters 
  theta <- runif(n, 0, 2*pi)
  phi <- runif(n, 0, pi/3)
  
  # noisy data 
  s <- 0.1
  x <- cos(theta) * sin(phi) + rnorm(n, 0, s)
  y <- sin(theta) * sin(phi) + rnorm(n, 0, s)
  z <- cos(phi) + rnorm(n, 0, s)
  
  return(list(x=x, y=y, z=z))
}

info <- function(D, method) {
  # model 
  m <- gam(z ~ te(x, y, bs = 'gp'), method=method, data=D)
  
  # details 
  pred <- predict(m, type='response', se.fit=TRUE)
  bias <- mean(pred$fit - D$z)
  variance <- mean(pred$se.fit^2)
  mse <- mean((D$z - pred$fit)^2)
  
  return(c(bias, variance, mse))
}

D <- surface.data(100)

GCV.Cp <- c(0, 0, 0) 
REML <- c(0, 0, 0)

# computational complexity, bias, variance, and mse
m <- microbenchmark(GCV.Cp <- GCV.Cp + info(D, method='GCV.Cp'), 
                    REML <-REML + info(D, method='REML'), 
                    times=250,
                    unit='ms')

# results 
m
GCV.Cp/250
REML/250
