# clear 
if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat("\014")

# set-up
library(splines)
set.seed(13)
heart <- read.table('SAheart.data.txt', sep=',', head=TRUE, row.names=1)

# sort data so it plots line left to right 
ind <-  sort(heart$tobacco, index.return=TRUE)$ix
heart$tobacco <- heart$tobacco[ind]
heart$chd <- heart$chd[ind]
x <- heart$tobacco

# bases
b.basis <- bs(x, df=6)
nat.basis <- ns(x, df=4)
trunc.basis.f <- function(x) {
  cbind(sapply(x, function(x) x), 
        sapply(x, function(x) x^2),
        sapply(x, function(x) x^3),
        sapply(x, function(x) pmax((x - 0.0525)^3, 0)),
        sapply(x, function(x) pmax((x - 2)^3, 0)),
        sapply(x, function(x) pmax((x - 5.5)^3, 0)))
}
trunc.basis <- trunc.basis.f(x)

# models 
b.model <- glm(chd ~ b.basis, data=heart, family=binomial)
nat.model <- glm(chd ~ nat.basis, data=heart, family=binomial)
trunc.model <- glm(chd ~ trunc.basis, data=heart, family=binomial)

# fitted values 
b.fit <- b.basis %*% coef(b.model)[2:7]
nat.fit <- nat.basis %*% coef(nat.model)[2:5]
trunc.fit <- trunc.basis %*% coef(trunc.model)[2:7]

# variances
b.H <- cbind(rep(1, 462), b.basis)
b.W <- diag(b.model$weight)
b.var <- solve(t(b.H) %*% b.W %*% b.H)[2:7, 2:7]
b.band <- sqrt(diag(b.basis %*% b.var %*% t(b.basis)))

nat.H <- cbind(rep(1, 462), nat.basis)
nat.W <- diag(nat.model$weight)
nat.var <- solve(t(nat.H) %*% nat.W %*% nat.H)[2:5, 2:5]
nat.band <- sqrt(diag(nat.basis %*% nat.var %*% t(nat.basis)))

trunc.H <- cbind(rep(1, 462), trunc.basis)
trunc.W <- diag(trunc.model$weight)
trunc.var <- solve((t(trunc.H) %*% trunc.W %*% trunc.H)+1e-6*diag(7))[2:7, 2:7]
trunc.band <- sqrt(diag(trunc.basis %*% trunc.var %*% t(trunc.basis)))

# plots
par(mfrow=c(3,1))

# B-Spline 
plot(x, b.fit, col='green', type='l', main='B-Spline', xlab='Tobacco', 
     ylab='Prediction', ylim=c(-1, 10), cex.lab=1.5)
polygon(c(x, rev(x)), c(b.fit - b.band, rev(b.fit + b.band)), col='yellow')
lines(x, b.fit, col='green', type='l')
axis(side=1, tck=0.05, at=x, labels=FALSE, tick=TRUE, col.ticks='red')

# Natural Spline 
plot(x, nat.fit, col='green', type='l', main='Natural Spline', xlab='Tobacco', 
     ylab='Prediction', ylim=c(-1, 10), cex.lab=1.5)
polygon(c(x, rev(x)), c(nat.fit-nat.band, rev(nat.fit+nat.band)), col='yellow')
lines(x, nat.fit, col='green', type='l')
axis(side=1, tck=0.05, at=x, labels=FALSE, tick=TRUE, col.ticks='red')

# Truncated Polynomial Spline 
plot(x, trunc.fit, col='green', type='l', main='Truncated Polynomial Spline', 
     xlab='Tobacco', ylab='Prediction', ylim=c(-1, 10), cex.lab=1.5)
polygon(c(x, rev(x)), c(trunc.fit - trunc.band, rev(trunc.fit + trunc.band)), 
        col='yellow')
lines(x, trunc.fit, col='green', type='l')
axis(side=1, tck=0.05, at=x, labels=FALSE, tick=TRUE, col.ticks='red')
