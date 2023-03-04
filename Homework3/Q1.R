# clear 
if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

# seed
set.seed(13)

# data for true curve - let's say sin(3x)/2
x.true <- seq(0.25, 2.6, length.out=100)
y.true <- sin(3*x.true)/2

# knots and noise
xi.1 <- 0.95
xi.2 <- 1.9
s <- 0.2

# generate as many points as the book did for each section
n1 <- 20
n2 <- 16
n3 <- 14

x1 <- runif(n1, min=0.25, max=xi.1) 
y1 <- sin(3*x1)/2 + rnorm(n1, mean=0, sd=s)

x2 <- runif(n2, min=xi.1, max=xi.2) 
y2 <- sin(3*x2)/2 + rnorm(n2, mean=0, sd=s)

x3 <- runif(n3, min=xi.2, max=2.6)
y3 <- sin(3*x3)/2 + rnorm(n3, mean=0, sd=s)

ys <- c(y1, y2, y3)
xs <- c(x1, x2, x3)


# plots the sample points and the true curve for the first three subplots   
plot.base <- function(t) {
  # figure set-up
  plot(x.true, y.true, xaxt='n', xlab='', ylab='', yaxt='n', main=t, 
       cex.main=1.5, ylim=c(min(ys) - 0.025, max(ys) + 0.025),
       xlim=c(min(xs) - 0.025, max(xs) + 0.025))
  
  rect(0.25, min(ys), 2.6, max(ys), col='#FDEB98', border=NA)
  rect(xi.1, min(ys), xi.2, max(ys), col='#CCBD7A', border=NA)
  
  # true curve and sample points 
  lines(x.true, y.true, col='blue', lwd=2) 
  points(xs, ys, lwd=1.25) 
  
  # dashed lines at knots and axis 
  abline(v=xi.1, col='black', lty=5) 
  abline(v=xi.2, col='black', lty=5)
  axis(1, at=c(xi.1, xi.2), labels=c(expression(xi[1]), expression(xi[2])), 
       cex.axis=1.5)
}

# basis functions according to ESL 
piecewise.constant <- function(x) {
  cbind(sapply(x, function(x) ifelse(x < xi.1, 1, 0)), 
        sapply(x, function(x) ifelse(x >= xi.1 & x < xi.2 , 1, 0)),
        sapply(x, function(x) ifelse(x >= xi.2, 1, 0)))
}

piecewise.linear <- function(x) {
  cbind(sapply(x, function(x) ifelse(x < xi.1, 1, 0)), 
        sapply(x, function(x) ifelse(x < xi.1, 1, 0)*x),
        sapply(x, function(x) ifelse(x >= xi.1 & x < xi.2 , 1, 0)),
        sapply(x, function(x) ifelse(x >= xi.1 & x < xi.2 , 1, 0)*x),
        sapply(x, function(x) ifelse(x >= xi.2, 1, 0)),
        sapply(x, function(x) ifelse(x >= xi.2, 1, 0)*x))
}

continuous.piecewise.linear <- function(x) {
  cbind(sapply(x, function(x) 1), 
        sapply(x, function(x) x),
        sapply(x, function(x) pmax(x-xi.1, 0)),
        sapply(x, function(x) pmax(x-xi.2, 0)))
}

spline.model <- function(B) {
  model <- lm(y~.-1, data=data.frame(x=B, y=ys))
  return(model$coefficients)
}


# 2 by 2 figure 
par(mfrow=c(2,2))


# plot (1, 1)
plot.base('Piecewise Constant')

# model 
beta.pc <- spline.model(piecewise.constant(xs))

# lines 
segments(x0=0.25, x1=xi.1, y0=beta.pc[1], y1=beta.pc[1], col='green', lwd=2)
segments(x0=xi.1, x1=xi.2, y0=beta.pc[2], y1=beta.pc[2], col='green', lwd=2)
segments(x0=xi.2, x1=2.6, y0=beta.pc[3], y1=beta.pc[3], col='green', lwd=2)


# plot (1, 2)
plot.base('Piecewise Linear')

# model 
beta.pl <- spline.model(piecewise.linear(xs))

# lines
segments(x0=0.25, x1=xi.1, 
         y0=beta.pl[1] + 0.25*beta.pl[2], 
         y1=beta.pl[1] + xi.1*beta.pl[2], 
         col='green', lwd=2)

segments(x0=xi.1, x1=xi.2, 
         y0=beta.pl[3] + xi.1*beta.pl[4], 
         y1=beta.pl[3] + xi.2*beta.pl[4], 
         col='green', lwd=2)

segments(x0=xi.2, x1=2.6, 
         y0=beta.pl[5] + xi.2*beta.pl[6], 
         y1=beta.pl[5] + 2.6*beta.pl[6], 
         col='green', lwd=2)


# plot (2, 1)
plot.base('Continuous Piecewise Linear')

# model 
beta.cpl <- spline.model(continuous.piecewise.linear(xs))

# lines
f1 <- beta.cpl[1] + beta.cpl[2]*0.25
f2 <- beta.cpl[1] + beta.cpl[2]*xi.1
f3 <- beta.cpl[1] + beta.cpl[2]*xi.2 + beta.cpl[3]*(xi.2 - xi.1)
f4 <- beta.cpl[1] + beta.cpl[2]*2.6 + beta.cpl[3]*(2.6 - xi.1) + 
  beta.cpl[4]*(2.6 - xi.2)

segments(x0=0.25, x1=xi.1, y0=f1, y1=f2, col='green', lwd=2)
segments(x0=xi.1, x1=xi.2, y0=f2, y1=f3, col='green', lwd=2)
segments(x0=xi.2, x1=2.6, y0=f3, y1=f4, col='green', lwd=2)


# plot (2, 2)

# basis function 
ys.b <- pmax(xs - xi.1, 0)

# plot base 
plot(xs, ys.b, xaxt='n', xlab='', ylab='', yaxt='n',
     main='Piecewise-linear Basis Function', col='white', cex.main=1.5, 
     ylim=c(min(ys.b) - 0.075, max(ys.b) + 0.075),
     xlim=c(min(xs) - 0.025, max(xs) + 0.025))
rect(0.25, 0, 2.6,  1.65, col='#FDEB98', border=NA)
rect(xi.1, 0, xi.2,  1.65, col='#CCBD7A', border=NA)
abline(v=xi.1, col='black', lty=5) 
abline(v=xi.2, col='black', lty=5)
text(1.4, 1, labels=expression((X- xi[1])['+']), cex=1.5)
axis(1, at=c(xi.1, xi.2), labels=c(expression(xi[1]), expression(xi[2])), 
     cex.axis=1.5)

# lines and points  
segments(x0=0.25, x1=xi.1, y0=0, y1=0, col='green', lwd=2)
segments(x0=xi.1, x1=2.6, y0=0, y1=pmax(2.6 - xi.1, 0), col='green', lwd=2)
points(xs, ys.b, cex=0.7, pch=16)
