if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat('\014')
set.seed(13)
library(ggplot2)
library(microbenchmark)
library(peakRAM)

nipals <- function(X, pcs = min(nrow(X), ncol(X)), 
                   center = TRUE, scale = TRUE, 
                   tol = 1e-6, max.iter = 500) {
  # ARGUMENTS 
  # X        - data (matrix)
  # pcs      - number of principal components (int)
  # center   - center the data to mean zero (bool)
  # scale    - scale data to unit variance (bool)
  # tol      - tolerance for convergence (numeric)
  # max.iter - max iterations without convergence (int)
  
  # VALUES   
  # P         - loadings (matrix)
  # Tm        - scores (matrix)
  # eig       - eigenvalues (vector)
  # explained - variance explained by each component (vector)

  # initialize        
  n <- nrow(X)                                   # rows
  p <- ncol(X)                                   # columns
  X <- scale(X, center = center, scale = scale)  # standardize 
  Tm <- matrix(0, nrow = n, ncol = pcs)          # score matrix 
  P <- matrix(0, nrow = p, ncol = pcs)           # loading matrix
  R <- X                                         # residual matrix
  TSS <- sum(R^2)                                # total sum of squares 
  eig <- rep(0, pcs)                             # eigenvalues 
  explained <- rep(0, pcs)                       # variance explained 
  
  # for each component 
  for (i in 1:pcs) {
    t.h <- R[, i] # initial score vector  

    # iterations 
    for (j in 1:max.iter) {
      # loading vector / GS correction / normalize 
      p.h <- crossprod(R, t.h) / sum(t.h^2)  
      if (i > 1){                                
        p.h <- p.h - P[, 1:(i-1)] %*% crossprod(P[, 1:(i-1)], p.h) 
      }
      p.h <- p.h / sqrt(sum(p.h^2))  
      
      # score vector / GS correction / score vector difference  
      t.h2 <- R %*% p.h   
      if (i > 1){
        t.h2 <- t.h2 - Tm[, 1:(i-1)] %*% crossprod(Tm[, 1:(i-1)], t.h2) 
      }
      e <- t.h2 - t.h
      
      # convergence check 
      if (sum(e^2) < tol){
        break
      } else{
        t.h <- t.h2
      }
    }
    
    # update and deflate 
    Tm[, i] <- t.h
    P[, i] <- p.h
    eig[i] <- sum(t.h^2)
    explained <- eig[i] / TSS
    R <- R - tcrossprod(t.h, p.h)
  }
  return(list(P=P, Tm=Tm, eig=eig, explained=explained))
}

# generates (n x p) data matrix
generate.data <- function(n, p){
  X <- scale(matrix(rnorm(n*p), nrow = n, ncol = p), center=TRUE, scale=TRUE)
  return(X)
}

# test sizes 
sizes <- c(10, 50, 100, 500, 1000, 1500, 2000)
N <- length(sizes)

# time benchmark 
time.list <- vector('list', N)
for (i in 1:N) {
    X <- generate.data(sizes[i], sizes[i])
    times <- microbenchmark(
      prcomp(X, center=FALSE, scale=FALSE),
      nipals(X, 2, center=FALSE, scale=FALSE),
      times = 10,
      unit ='ms')
    # store 
    time.list[[i]] <- summary(times)$mean
    print(summary(times)$mean)
}

# time data  
time.r <- log10(sapply(time.list, '[', 1))
time.nipals <- log10(sapply(time.list, '[', 2))
logx <- log10(sizes)

# time plot 
p <- ggplot() +
  geom_line(aes(x=logx, y=time.r, color='SVD'), size = 1.15) +
  geom_point(aes(x=logx, y=time.r), color='black') +
  geom_line(aes(x=logx, y=time.nipals, color='NIPALS'), size = 1.15) + 
  geom_point(aes(x=logx, y=time.nipals), color='black') +
  scale_color_manual(values = c('SVD'='cyan', 'NIPALS'='darkorange')) +
  theme(text = element_text(size=20)) +
  labs(x = 'log-size', 
       y = 'log-time (ms)',
       title = 'SVD vs NIPALS runtime',
       color = 'Method')  

# data
X <- generate.data(10^3, 10^3)

# accuracy and memory
peak.r <- 0
peak.nipals <- 0
accuracy <- 0 
for (i in 1:10) {
  mem.r <- peakRAM(P.r <- prcomp(X, scale = FALSE, center = FALSE))
  mem.nipals <- peakRAM(P.nipals <- nipals(X, pcs = 2, 
                                           scale = FALSE, center = FALSE))
  if (mem.r$Peak_RAM_Used_MiB > peak.r) {
    peak.r <- mem.r$Peak_RAM_Used_MiB
  }
  if (mem.nipals$Peak_RAM_Used_MiB > peak.nipals) {
    peak.nipals <- mem.nipals$Peak_RAM_Used_MiB
  accuracy <- accuracy + norm(abs(P.r$rotation[, 1:2]) - abs(P.nipals$P), "2")
  }
}
accuracy <- accuracy/10

peak.r
peak.nipals
accuracy
p