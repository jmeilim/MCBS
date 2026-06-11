library(Ball)




MCBS.cor <- function(x, y, A, M, distance = FALSE) {
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  idx00 <- which(A == 0 & M == 0)
  idx01 <- which(A == 0 & M == 1)
  idx10 <- which(A == 1 & M == 0)
  idx11 <- which(A == 1 & M == 1)
  
  p00 <- length(idx00)/length(A)
  p01 <- length(idx01)/length(A)
  p10 <- length(idx10)/length(A)
  p11 <- length(idx11)/length(A)
  
  if(distance == TRUE){
    
    x00 <- x[idx00, idx00]
    y00 <- y[idx00, idx00]
    
    x01 <- x[idx01, idx01]
    y01 <- y[idx01, idx01]
    
    x10 <- x[idx10, idx10]
    y10 <- y[idx10, idx10]
    
    x11 <- x[idx11, idx11]
    y11 <- y[idx11, idx11]
    
    MCBS.cor <-
      p00*bcov(x00,y00,distance=TRUE)^2 +
      p01*bcov(x01,y01,distance=TRUE)^2 +
      p10*bcov(x10,y10,distance=TRUE)^2 +
      p11*bcov(x11,y11,distance=TRUE)^2
    
  } else {
    
    x00 <- x[idx00, ]
    y00 <- y[idx00, ]
    
    x01 <- x[idx01, ]
    y01 <- y[idx01, ]
    
    x10 <- x[idx10, ]
    y10 <- y[idx10, ]
    
    x11 <- x[idx11, ]
    y11 <- y[idx11, ]
    
    MCBS.cor <-
      p00*bcov(x00,y00)^2 +
      p01*bcov(x01,y01)^2 +
      p10*bcov(x10,y10)^2 +
      p11*bcov(x11,y11)^2
  }
  
  return(MCBS.cor)
}

Causal.cor <- function(x, y, z, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha = length(index0)/length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]
    y0 <- y[index0, index0]
    x1 <- x[index1, index1]
    y1 <- y[index1, index1]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, x0, distance = TRUE)^2 + (1-alpha)*bcov(x1, y1, distance = TRUE)^2
  } else {
    x0 <- x[index0, ]
    y0 <- y[index0, ]
    x1 <- x[index1, ]
    y1 <- y[index1, ]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, y0)^2 + (1-alpha)*bcov(x1, y1)^2
  }
  return(Causal.cor)
}


screening <- function(X, Y, M, A, k) {
  p <- ncol(X)
  S1 <- S2 <- S3 <- numeric(p)
  
  for(j in 1:p) {
    xj <- X[, j]
    S1[j] <- Causal.cor(xj, Y, A)
    S2[j] <- Causal.cor(xj, M, A)
    S3[j] <- MCBS.cor(as.matrix(xj), as.matrix(Y), A, M)
  }
  
  r1 <- rank(-S1)
  r2 <- rank(-S2)
  r3 <- rank(-S3)
  geo <- (r1 * r2 * r3)^(1/3)
  
  res <- data.frame(
    variable = 1:p,
    S1 = S1,
    S2 = S2,
    S3 = S3,
    r1 = r1,
    r2 = r2,
    r3 = r3,
    geo = geo
  )
  
  res <- res[order(res$geo), ]
  
  print(res[1:10, ])
  cat("Rangs des vrais confondants :\n")
  print(res[res$variable %in% c(1,2,3), ])
  
  return(order(geo)[1:k])
}


S1_score_clean <- function(xj, Y, A) {
  w1 <- mean(A == 1)
  w0 <- mean(A == 0)
  
  s1 <- 0
  
  if (sum(A == 1) >= 2) {
    s1 <- s1 + w1 * bcov(xj[A == 1], Y[A == 1])^2
  }
  
  if (sum(A == 0) >= 2) {
    s1 <- s1 + w0 * bcov(xj[A == 0], Y[A == 0])^2
  }
  
  return(s1)
}

screening_S1_clean <- function(X, Y, A, k) {
  p <- ncol(X)
  S1 <- numeric(p)
  
  for (j in 1:p) {
    S1[j] <- S1_score_clean(X[, j], Y, A)
  }
  
  return(order(-S1)[1:k])
}

