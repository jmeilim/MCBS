# data generation process


expit <- function(x) {
  exp(x) / (1 + exp(x))
}

DGP <- function(n, p, beta_AY,  beta_MY, beta_AM) { 
  
  X <- matrix(runif(n*p, -1, 1), ncol = p)
  
  
  prob_A <- expit(X[,1]*beta_AY +
                    X[,2]*beta_AM)
  A <- rbinom(n, 1, prob_A)
  
  
  prob_M <- expit(A +
                    X[,2]*beta_AM +
                    X[,3]*beta_MY)
  M <- rbinom(n, 1, prob_M)
  
  
  Y <- 1.5*A + 2*M +
    X[,1]*beta_AY +
    X[,3]*beta_MY +
    rnorm(n)
  
  return(list(X=X, A=A, M=M, Y=Y))
}

DGP_rich <- function(
    n,
    p,
    lambda = 1,
    beta_base = c(weak = 0.1, medium = 0.2, strong = 0.3)
) {
  
  X <- matrix(runif(n * p, -1, 1), ncol = p)
  
  beta <- lambda * beta_base
  
 
  idx_AY <- 1:3
  idx_AM <- 4:6
  idx_MY <- 7:9
  idx_pY <- 10:12
  idx_pM <- 13:15
  idx_pA <- 16:20
  
  b_AY <- beta
  b_AM <- beta
  b_MY <- beta
  b_pY <- beta
  b_pM <- beta
  b_pA <- c(0.1 , 0.2 , 0.4 , 0.6 , 0.8)
  
  eta_A <-
    X[, idx_AY, drop = FALSE] %*% b_AY +
    X[, idx_AM, drop = FALSE] %*% b_AM +
    X[, idx_pA, drop = FALSE] %*% b_pA
  
  prob_A <- expit(as.vector(eta_A))
  A <- rbinom(n, 1, prob_A)
  
  eta_M <-
    A +
    X[, idx_AM, drop = FALSE] %*% b_AM +
    X[, idx_MY, drop = FALSE] %*% b_MY +
    X[, idx_pM, drop = FALSE] %*% b_pM
  
  prob_M <- expit(as.vector(eta_M))
  M <- rbinom(n, 1, prob_M)
  
  Y <-
    1.5 * A +
    2 * M +
    X[, idx_AY, drop = FALSE] %*% b_AY +
    X[, idx_MY, drop = FALSE] %*% b_MY +
    X[, idx_pY, drop = FALSE] %*% b_pY +
    rnorm(n)
  
  Y <- as.vector(Y)
  
  true_confounders <- c(idx_AY, idx_AM, idx_MY)
  pure_predictors <- c(idx_pY, idx_pM)
  true_signals <- c(true_confounders, pure_predictors)
  
  variable_info <- data.frame(
    variable = true_signals,
    type = c(
      rep("AY", 3),
      rep("AM", 3),
      rep("MY", 3),
      rep("pY", 3),
      rep("pM", 3)
    ),
    level = rep(c("weak", "medium", "strong"), times = 5),
    beta = rep(as.numeric(beta), times = 5)
  )
  
  return(list(
    X = X,
    A = A,
    M = M,
    Y = Y,
    true_confounders = true_confounders,
    pure_predictors = pure_predictors,
    true_signals = true_signals,
    variable_info = variable_info,
    beta = beta
  ))
}

