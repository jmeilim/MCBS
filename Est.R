
library(parallel)
library(Ball)
library(glmnet)


N       <- 300      
P       <- 1000     
ALPHA   <- 1.0      
BETA_C  <- 2        
BETA_D  <- 10       
TRUTH   <- 2        
B_PERM  <- 100      
GAMMA   <- 0.10    
N_RUNS  <- 500     
NCORES  <- 120
TRIM    <- 0.02
CONF    <- c(1, 2)      
TARGET  <- c(1, 2, 4)     




generer <- function(seed) {
  set.seed(seed)
  X  <- matrix(runif(N * P, -1, 1), ncol = P)
  ps <- 1 / (1 + exp(-rowSums(X[, CONF]) * ALPHA))
  D  <- rbinom(N, 1, ps)
  Y  <- D * TRUTH + rowSums(X[, CONF]) * BETA_C + X[, 4] * BETA_D + rnorm(N)
  list(X = X, D = D, Y = Y)
}



score_propension <- function(X, D, K) {
  Xs <- as.data.frame(X[, K, drop = FALSE])
  names(Xs) <- paste0("v", seq_along(K))
  fit <- suppressWarnings(
    glm(D ~ ., data = data.frame(D = D, Xs), family = binomial))
  ps <- as.numeric(predict(fit, type = "response"))
  pmin(pmax(ps, TRIM), 1 - TRIM)
}



estim_ipw <- function(X, D, Y, K, hajek = TRUE) {
  if (length(K) == 0) return(NA_real_)
  ps <- score_propension(X, D, K)
  w  <- D / ps + (1 - D) / (1 - ps)
  if (hajek) {
    m1 <- sum(w[D == 1] * Y[D == 1]) / sum(w[D == 1])
    m0 <- sum(w[D == 0] * Y[D == 0]) / sum(w[D == 0])
  } else {
    m1 <- sum(w[D == 1] * Y[D == 1]) / length(Y)
    m0 <- sum(w[D == 0] * Y[D == 0]) / length(Y)
  }
  m1 - m0
}



estim_dr <- function(X, D, Y, K) {
  if (length(K) == 0) return(NA_real_)
  Xs <- as.data.frame(X[, K, drop = FALSE])
  names(Xs) <- paste0("v", seq_along(K))
  ps <- score_propension(X, D, K)
  
  predire_bras <- function(idx) {
    m <- lm(Y ~ ., data = data.frame(Y = Y[idx], Xs[idx, , drop = FALSE]))
    as.numeric(predict(m, newdata = Xs))
  }
  b1 <- predire_bras(which(D == 1))
  b0 <- predire_bras(which(D == 0))
  
  mean( (D * Y - (D - ps) * b1) / ps -
          ((1 - D) * Y + (D - ps) * b0) / (1 - ps) )
}



replication <- function(seed) {
  d <- generer(seed)
  X <- d$X; D <- d$D; Y <- d$Y
  

  scores <- vapply(1:P, function(j) Causal.cor(X[, j], Y, D), numeric(1))
  ordre  <- order(scores, decreasing = TRUE)
  K_cbs  <- ordre[1:30]
  
  
  peel   <- peel_until_noise(X, D, Y, gamma = GAMMA, B = B_PERM, verbose = FALSE)
  K_peel <- peel$selected
  
  c(dr_oracle          = estim_dr (X, D, Y, CONF),
    dr_oracle_plus     = estim_dr (X, D, Y, TARGET),
    dr_cbs             = estim_dr (X, D, Y, K_cbs),
    dr_peel            = estim_dr (X, D, Y, K_peel),
    
    ipwHT_oracle       = estim_ipw(X, D, Y, CONF,   hajek = FALSE),
    ipwHT_oracle_plus  = estim_ipw(X, D, Y, TARGET, hajek = FALSE),
    ipwHT_cbs          = estim_ipw(X, D, Y, K_cbs,  hajek = FALSE),
    ipwHT_peel         = estim_ipw(X, D, Y, K_peel, hajek = FALSE),
    
    ipwHJ_oracle       = estim_ipw(X, D, Y, CONF,   hajek = TRUE),
    ipwHJ_oracle_plus  = estim_ipw(X, D, Y, TARGET, hajek = TRUE),
    ipwHJ_cbs          = estim_ipw(X, D, Y, K_cbs,  hajek = TRUE),
    ipwHJ_peel         = estim_ipw(X, D, Y, K_peel, hajek = TRUE),
    
    cap_cbs      = as.numeric(all(CONF %in% K_cbs)),
    cap_peel     = as.numeric(all(CONF %in% K_peel)),
    cible_exacte = as.numeric(setequal(K_peel, TARGET)),
    k_hat        = peel$k_hat,
    rang_X1      = match(1, ordre),
    rang_X2      = match(2, ordre),
    rang_X4      = match(4, ordre))
}



cat(sprintf("Lancement : %d runs, n=%d, p=%d, alpha=%.1f, beta_dom=%g, B=%d\n\n",
            N_RUNS, N, P, ALPHA, BETA_D, B_PERM))

cl <- makeCluster(NCORES)
clusterEvalQ(cl, { library(Ball); library(glmnet) })
clusterExport(cl, c("generer", "score_propension", "estim_ipw", "estim_dr",
                    "replication", "peel_until_noise", "noise_threshold",
                    "permute_within_A", "resid_on", "Causal.cor", "expit",
                    "N", "P", "ALPHA", "BETA_C", "BETA_D", "TRUTH",
                    "B_PERM", "GAMMA", "TRIM", "CONF", "TARGET"),
              envir = .GlobalEnv)

duree <- system.time(
  MC <- as.data.frame(t(parSapply(cl, 1:N_RUNS, replication)))
)["elapsed"]
stopCluster(cl)

cat(sprintf("Termine en %.1f min\n\n", duree / 60))



resumer <- function(v) c(
  biais = mean(v) - TRUTH,
  mc_se = sd(v) / sqrt(length(v)),
  sd    = sd(v),
  rmse  = sqrt(mean((v - TRUTH)^2))
)

estimateurs <- c("dr_oracle",    "dr_oracle_plus",    "dr_cbs",    "dr_peel",
                 "ipwHT_oracle", "ipwHT_oracle_plus", "ipwHT_cbs", "ipwHT_peel",
                 "ipwHJ_oracle", "ipwHJ_oracle_plus", "ipwHJ_cbs", "ipwHJ_peel")
table_res <- t(sapply(estimateurs, function(nm) resumer(MC[[nm]])))

cat("===== ESTIMATION =====\n")
print(round(table_res, 4))

cat("\n===== SCREENING =====\n")
cat(sprintf("  capture {1,2}        : CBS %3.0f%%   |   peeling %3.0f%%\n",
            100 * mean(MC$cap_cbs), 100 * mean(MC$cap_peel)))
cat(sprintf("  cible exacte {1,2,4} : %3.0f%%   |   k_hat median %.0f [%.0f, %.0f]\n",
            100 * mean(MC$cible_exacte),
            median(MC$k_hat), min(MC$k_hat), max(MC$k_hat)))
cat(sprintf("  rangs BCov (medianes) : X1 %4.0f   X2 %4.0f   X4 %4.0f\n",
            median(MC$rang_X1), median(MC$rang_X2), median(MC$rang_X4)))

cat("\n===== HT vs HAJEK : gain de variance =====\n")
for (ens in c("oracle", "oracle_plus", "cbs", "peel")) {
  ht <- MC[[paste0("ipwHT_", ens)]]
  hj <- MC[[paste0("ipwHJ_", ens)]]
  cat(sprintf("  %-12s : sd  HT %.4f -> HJ %.4f  (%+.0f%%)   |  biais HT %+.4f  HJ %+.4f\n",
              ens, sd(ht), sd(hj), 100*(sd(hj)/sd(ht) - 1),
              mean(ht) - TRUTH, mean(hj) - TRUTH))
}
