library(energy) 


Causal.dcor <- function(x, y, A) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  i0 <- which(A == 0)
  i1 <- which(A == 1)
  alpha <- length(i0) / length(A)
  

  d0 <- dcor(x[i0, , drop = FALSE], y[i0, , drop = FALSE])^2
  d1 <- dcor(x[i1, , drop = FALSE], y[i1, , drop = FALSE])^2
  
  alpha * d0 + (1 - alpha) * d1
}

resid_Y <- function(Y, A, X, K) {
  as.numeric(residuals(lm(Y ~ A + X[, K, drop = FALSE])))
}

B <- 30
R <- t(sapply(1:B, function(s) {
  set.seed(s)
  X  <- matrix(runif(n*p, -1, 1), ncol = p)
  ps <- 1/(1 + exp(-rowSums(X[, 1:2])))
  A  <- rbinom(n, 1, ps)
  Y  <- A*2 + rowSums(X[, 1:2])*2 + X[, 4]*10 + rnorm(n)
  
 
  sb0 <- vapply(1:p, function(j) Causal.cor (X[, j], Y, A), numeric(1))
  sd0 <- vapply(1:p, function(j) Causal.dcor(X[, j], Y, A), numeric(1))
  

  Yr  <- resid_Y(Y, A, X, 4)
  sb1 <- vapply(1:p, function(j) Causal.cor (X[, j], Yr, A), numeric(1))
  sd1 <- vapply(1:p, function(j) Causal.dcor(X[, j], Yr, A), numeric(1))
  
  rg <- function(s, v) match(v, order(s, decreasing = TRUE))
  c(bcov_av_X1 = rg(sb0,1), bcov_av_X2 = rg(sb0,2),
    bcov_ap_X1 = rg(sb1,1), bcov_ap_X2 = rg(sb1,2),
    dcor_av_X1 = rg(sd0,1), dcor_av_X2 = rg(sd0,2),
    dcor_ap_X1 = rg(sd1,1), dcor_ap_X2 = rg(sd1,2))
}))

med <- function(col) median(R[, col])
cat("                 X1        X2\n")
cat(sprintf("BCov  avant :  %4.0f      %4.0f\n", med("bcov_av_X1"), med("bcov_av_X2")))
cat(sprintf("BCov  apres :  %4.0f      %4.0f\n", med("bcov_ap_X1"), med("bcov_ap_X2")))
cat(sprintf("DC    avant :  %4.0f      %4.0f\n", med("dcor_av_X1"), med("dcor_av_X2")))
cat(sprintf("DC    apres :  %4.0f      %4.0f\n", med("dcor_ap_X1"), med("dcor_ap_X2")))