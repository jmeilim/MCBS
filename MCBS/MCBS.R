# ============================================================
# MCBS.R — CBS Mediation
# ============================================================
library(Ball)
library(glmnet)


MCBS.cor <- function(x, y, A, M, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  idx00 <- which(A==0 & M==0); idx01 <- which(A==0 & M==1)
  idx10 <- which(A==1 & M==0); idx11 <- which(A==1 & M==1)
  p00 <- length(idx00)/length(A); p01 <- length(idx01)/length(A)
  p10 <- length(idx10)/length(A); p11 <- length(idx11)/length(A)
  if (distance == TRUE) {
    MCBS.cor <-
      p00*bcov(x[idx00,idx00], y[idx00,idx00], distance=TRUE)^2 +
      p01*bcov(x[idx01,idx01], y[idx01,idx01], distance=TRUE)^2 +
      p10*bcov(x[idx10,idx10], y[idx10,idx10], distance=TRUE)^2 +
      p11*bcov(x[idx11,idx11], y[idx11,idx11], distance=TRUE)^2
  } else {
    MCBS.cor <-
      p00*bcov(x[idx00,], y[idx00,])^2 +
      p01*bcov(x[idx01,], y[idx01,])^2 +
      p10*bcov(x[idx10,], y[idx10,])^2 +
      p11*bcov(x[idx11,], y[idx11,])^2
  }
  return(MCBS.cor)
}


Causal.cor <- function(x, y, z, distance = FALSE) {
  x      <- as.matrix(x)
  y      <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha  <- length(index0) / length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]; y0 <- y[index0, index0]
    x1 <- x[index1, index1]; y1 <- y[index1, index1]
    Causal.cor <- alpha*(1-alpha)*bcov(x0,y0,distance=TRUE)^2 +
      (1-alpha)*alpha*bcov(x1,y1,distance=TRUE)^2
  } else {
    x0 <- x[index0, ]; y0 <- y[index0, ]
    x1 <- x[index1, ]; y1 <- y[index1, ]
    Causal.cor <- alpha*bcov(x0,y0)^2 + (1-alpha)*bcov(x1,y1)^2
  }
  return(Causal.cor)
}


