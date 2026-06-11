# data generation process


expit <- function(x){
  exp(x)/(1+exp(x))
}

DGP <- function(n, p, beta_AY,  beta_MY, beta_AM) { 
  
  X <- matrix(runif(n*p, -1, 1), ncol = p)
  

  prob_A <- expit(X[,1]*beta_AY +
                  X[,2]*beta_AM + X[,4])
  A <- rbinom(n, 1, prob_A)
  

  prob_M <- expit(A +
                  X[,2]*beta_AM +
                  X[,3]*beta_MY)
  M <- rbinom(n, 1, prob_M)
  

  Y <- 1.5*A + 2*M +
       X[,1]*beta_AY +
       X[,3]*beta_MY +
       X[,5] +
       rnorm(n)
  
  return(list(X=X, A=A, M=M, Y=Y))
}

