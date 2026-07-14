# data generation process
library(MASS)

expit <- function(x) {
  exp(x) / (1 + exp(x))
}




expit <- function(x) 1 / (1 + exp(-x))

DGP_oal_total <- function(n, p, rho_x = 0, sig_e = 0.001, bA = 1,
                          beta_strong = 1.5, beta_weak = 1,
                          alpha_conf = 1.0, alpha_instr = 1.0) {
  
  # -- 6 covariees non-bruit, 2 par role, comme OAL
  pC <- 2  # confondeurs FORTS   : Xc (beta fort, alpha fort)
  pP <- 2  # confondeurs FAIBLES : Xp (beta faible, alpha fort)  <- domination
  pI <- 2  # instruments         : Xi (beta = 0,   alpha fort)
  pS <- p - (pC + pP + pI)
  
  # -- coefficients
  beta_v  <- c(rep(beta_strong, pC),
               rep(beta_weak,   pP),
               rep(0,           pI),
               rep(0,           pS))
  
  alpha_v <- c(rep(alpha_conf,  pC),
               rep(alpha_conf,  pP),
               rep(alpha_instr, pI),
               rep(0,           pS))
  
  # -- X ~ N(0, Sigma) compound symmetry (rho_x hors diagonale)
  Sigma <- matrix(rho_x, nrow = p, ncol = p)
  diag(Sigma) <- 1
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # -- traitement
  pA <- expit(as.vector(X %*% alpha_v))
  A  <- rbinom(n, 1, pA)
  
  # -- outcome (effet total seul, PAS de M)
  Y <- as.vector(X %*% beta_v) + bA * A + rnorm(n, sd = sig_e)
  
  # -- standardisation post-generation (comme OAL)
  X <- scale(X)
  
  # -- meta compatible avec tes autres DGP
  AY_strong <- 1:pC                          # Xc : confondeurs forts
  AY_weak   <- (pC + 1):(pC + pP)            # Xp : confondeurs faibles
  Instr     <- (pC + pP + 1):(pC + pP + pI)  # Xi : instruments
  
  true_confounders <- c(AY_strong, AY_weak)  # ce qui compte pour l'ATE
  true_signals     <- c(true_confounders, Instr)
  
  variable_info <- data.frame(
    variable = true_signals,
    type     = c(rep("AY", pC + pP), rep("Instr", pI)),
    level    = c(rep("strong", pC), rep("weak", pP), rep("strong", pI)),
    beta     = c(rep(beta_strong, pC), rep(beta_weak, pP), rep(0, pI)),
    alpha    = c(rep(alpha_conf, pC + pP), rep(alpha_instr, pI))
  )
  
  list(
    X = X, A = A, Y = Y,
    true_signals     = true_signals,
    true_confounders = true_confounders,
    variable_info    = variable_info
  )
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


  
DGP_four_signals_corr<- function(n, p, rho_x = 0, noise_sd = 0.3) {
    
  expit <- function(x) {
    1 / (1 + exp(-x))
  }
  
  # ============================================================
  # Coefficients définis une seule fois
  # ============================================================
  
  beta_AY <- c(
    weak   = 0.1,
    strong = 1.00
  )
  
  beta_MY <- c(
    weak   = 0.2,
    strong = 1.00
  )
  
  beta_A_to_Y <- 1.5
  beta_M_to_Y <- 2.0
  beta_A_to_M <- 1.0
  
  # ============================================================
  # Génération de X avec corrélation rho_x
  # ============================================================
  
  common_factor <- rnorm(n)
  E <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  X <- sqrt(rho_x) * common_factor + sqrt(1 - rho_x) * E
  
  # ============================================================
  # Variables actives
  # ============================================================
  
  AY_vars <- 1:2
  MY_vars <- 3:4
  
  AY <- X[, AY_vars, drop = FALSE]
  MY <- X[, MY_vars, drop = FALSE]
  
  colnames(AY) <- names(beta_AY)
  colnames(MY) <- names(beta_MY)
  
  # ============================================================
  # Traitement A : AY agit sur A
  # ============================================================
  
  eta_A <- AY %*% beta_AY
  
  prob_A <- expit(eta_A)
  A <- rbinom(n, 1, prob_A)
  
  # ============================================================
  # Médiateur M : MY agit sur M
  # ============================================================
  
  eta_M <- beta_A_to_M * A + MY %*% beta_MY
  
  prob_M <- expit(eta_M)
  M <- rbinom(n, 1, prob_M)
  
  # ============================================================
  # Outcome Y : AY et MY agissent sur Y
  # ============================================================
  
  eps <- rnorm(n, 0, noise_sd)
  
  Y <- beta_A_to_Y * A +
    beta_M_to_Y * M +
    AY %*% beta_AY +
    MY %*% beta_MY +
    eps
  
  Y <- as.numeric(Y)
  
  # ============================================================
  # Informations sur les vraies variables
  # ============================================================
  
  true_signals <- c(AY_vars, MY_vars)
  true_confounders <- true_signals
  
  variable_info <- data.frame(
    variable = true_signals,
    type = c(rep("AY", length(beta_AY)), rep("MY", length(beta_MY))),
    level = c(names(beta_AY), names(beta_MY)),
    beta = c(beta_AY, beta_MY)
  )
  
  list(
    X = X,
    A = A,
    M = M,
    Y = Y,
    true_signals = true_signals,
    true_confounders = true_confounders,
    variable_info = variable_info
  )
}


DGP_cbs_style <- function(n = 200, p = 1000,
                          bA = 2, beta_dom = 10, beta_weak = 1,
                          n_weak = 13, alpha_conf = 2, alpha_instr = 3) {
  
  X <- matrix(runif(n * p, min = -1, max = 1), ncol = p)
  
  # -- roles SANS chevauchement ------------------------------------
  idx_dom   <- 1                          # dominant sur Y (pas dans A)
  idx_conf  <- c(2, 3)                    # confondeurs (A et Y)
  # predicteurs purs de Y : juste apres les confondeurs, n_weak-2 d'entre eux
  idx_pureY <- 4:(1 + n_weak)             # X4..X14  (11 variables si n_weak=13)
  idx_instr <- (2 + n_weak):(3 + n_weak)  # X15, X16 : instruments, hors zone Y
  
  # variables agissant sur Y (dominant + confondeurs + purs)
  idx_Y     <- c(idx_dom, idx_conf, idx_pureY)
  
  # -- traitement : confondeurs + instruments agissent sur A -------
  eta_A <- rowSums(X[, idx_conf,  drop = FALSE]) * alpha_conf +
    rowSums(X[, idx_instr, drop = FALSE]) * alpha_instr
  A <- rbinom(n, 1, expit(eta_A))
  
  # -- outcome : dominant (beta_dom) + confondeurs & purs (beta_weak)
  Y <- bA * A +
    X[, idx_dom] * beta_dom +
    rowSums(X[, c(idx_conf, idx_pureY), drop = FALSE]) * beta_weak +
    rnorm(n)
  
  # -- meta (aucun chevauchement) ----------------------------------
  variable_info <- data.frame(
    variable = c(idx_dom, idx_conf, idx_pureY, idx_instr),
    type = c("domPred",
             rep("confounder", length(idx_conf)),
             rep("purePred",   length(idx_pureY)),
             rep("instrument", length(idx_instr))),
    level = c("dominant",
              rep("weak", length(idx_conf) + length(idx_pureY)),
              rep("none", length(idx_instr))),
    beta = c(beta_dom,
             rep(beta_weak, length(idx_conf) + length(idx_pureY)),
             rep(0, length(idx_instr))),
    stringsAsFactors = FALSE
  )
  
  list(
    X = X, A = A, Y = Y,
    true_confounders = idx_conf,
    pure_predictors  = idx_pureY,
    dominant         = idx_dom,
    instruments      = idx_instr,
    true_signals     = variable_info$variable,
    variable_info    = variable_info,
    bA = bA
  )
}