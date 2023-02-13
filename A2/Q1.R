if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

library(microbenchmark)
library(ggplot2)
library(tictoc)
set.seed(13)

# coefficients via naive linear algebra 
naive.lm <- function(X, Y){
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  return(coefficients=as.vector(beta))
}

# coefficients via QR decomposition
qr.lm <- function(X, Y){
  # X = QR 
  decomp <- qr(X)  
  
  # solve R*beta = Q'*Y
  beta <- solve.qr(decomp, Y)  
  return(coefficients=as.vector(beta)) 
}

# coefficients via SVD
svd.lm <- function(X, Y){
  # X = UDV'
  decomp <- svd(X)  
  U <- decomp$u  
  D.inv <- diag(1/decomp$d)
  V <- decomp$v 
    
  beta <- V %*% D.inv %*% t(U) %*% Y
  return(coefficients=as.vector(beta))
}

# generates (n x p) design matrix and (1 x n) response vector 
generate.data <- function(n, p){
  Y <- rnorm(n)
  X <- matrix(rnorm(p*n), ncol = p)
  list(X = X, Y = Y)
}

# set-up variables 
sizes <- c(100, 250, 500, 1000, 2500, 5000, 10^4, 10^5/4,
           10^5/2, 10^5, 10^6/4, 10^6/2, 10^6)
n <- length(sizes)
init <- rep(0, n)
lm.times <- init
naive.times <- init
qr.times <- init
svd.times <- init

# main loop
tic()
for (i in 1:n){
  # benchmark 
  D <- generate.data(sizes[i], p = 10)
  times <- microbenchmark(
    as.vector(lm.fit(D$X, D$Y)$coefficients),
    naive.lm(D$X, D$Y),
    qr.lm(D$X, D$Y), 
    svd.lm(D$X, D$Y),
    times= 50,
    check='equal',
    unit='ms')
  s <- summary(times)$mean
  
  # append
  lm.times[i] <- s[1]
  naive.times[i] <- s[2]
  qr.times[i] <- s[3]
  svd.times[i] <- s[4]
}

# plotting variables 
times.df <- list(data.frame(x=sizes, y=lm.times), 
                 data.frame(x=sizes, y=naive.times),
                 data.frame(x=sizes, y=qr.times),
                 data.frame(x=sizes, y=svd.times))
c <- c('black', 'red', 'blue', 'green')

plot.times <- function(data, c){
  # set-up ggplot
  pl <- ggplot() + 
    geom_point(data = times.df[[1]], aes(x, y, color = 'black')) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('Lienar Regression Coefficients Benchmarking') +
    xlab('log(n)') + 
    ylab('log(time) in ms') +
    scale_color_manual(name='Methods',
                       values=c('Base R'=c[1], 
                                'Naive Linear Algebra'=c[2], 
                                'QR Decomposition'=c[3],
                                'SVD'=c[4]))
  
  # plot lines and points
  for (i in 1:length(data)){
    pl <- pl +
      geom_line(data = data[[i]], aes(x, y, color = c[i]), color=c[i]) + 
      geom_point(data = data[[i]], aes(x, y, color = c[i]), color=c[i])
  }
  return(pl)
}

plot.times(times.df, c)

# log-log models 
model.lm <- lm(log(y) ~ log(x), data=times.df[[1]])
model.naive <- lm(log(y) ~ log(x), data=times.df[[2]])
model.qr <- lm(log(y) ~ log(x), data=times.df[[3]])
model.svd <- lm(log(y) ~ log(x), data=times.df[[4]])

model.lm
model.naive
model.qr
model.svd
toc()  # 158.244 seconds 