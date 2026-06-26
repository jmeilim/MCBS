# ============================================================
# 10 geo13 raw + 10 S2 raw + 10 geo13 residualized Y
# B = 100 Monte Carlo
# ============================================================

# Suppose déjà définis :
# DGP_rich
# compute_scores
# Causal.cor
# MCBS.cor

# ============================================================
# 1. Résidu simple de Y
# ============================================================

resid_lm <- function(y, Z) {
  Z <- as.data.frame(Z)
  df <- data.frame(y = y, Z)
  
  fit <- tryCatch(
    lm(y ~ ., data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(y - mean(y))
  }
  
  as.numeric(resid(fit))
}

# ============================================================
# 2. Scores après résidualisation de Y seulement
# ============================================================

compute_scores_resid_Yonly <- function(X, A, M, Y, K) {
  p <- ncol(X)
  
  Yres <- resid_lm(
    Y,
    cbind(A = A, M = M, X[, K, drop = FALSE])
  )
  
  S1 <- S2 <- S3 <- numeric(p)
  
  for (j in 1:p) {
    S1[j] <- Causal.cor(X[, j], Yres, A)
    S2[j] <- Causal.cor(X[, j], M, A)  # M brut, pas de résidualisation
    S3[j] <- MCBS.cor(as.matrix(X[, j]), as.matrix(Yres), A, M)
  }
  
  list(S1 = S1, S2 = S2, S3 = S3)
}

# ============================================================
# 3. Ranking geo13 = geo(S1,S3)
# ============================================================

order_geo13 <- function(sc, candidates) {
  r1 <- rank(-sc$S1[candidates], ties.method = "average")
  r3 <- rank(-sc$S3[candidates], ties.method = "average")
  
  score <- sqrt(r1 * r3)
  
  candidates[order(score)]
}

make_rank <- function(ord, p) {
  rk <- rep(NA, p)
  rk[ord] <- seq_along(ord)
  rk
}

# ============================================================
# 4. Paramètres Monte Carlo
# ============================================================

set.seed(42)

B <- 100
n <- 300
p <- 1000
lambda <- 4

k_raw_Y <- 10
k_raw_M <- 10
k_res_Y <- 10

res <- data.frame()

# ============================================================
# 5. Boucle Monte Carlo
# ============================================================

for (b in 1:B) {
  
  cat("Run", b, "/", B, "\n")
  
  out <- tryCatch({
    
    dat <- DGP_rich(n, p, lambda = lambda)
    
    X <- dat$X
    A <- dat$A
    M <- dat$M
    Y <- dat$Y
    
    true_active <- dat$true_signals
    true_conf <- dat$true_confounders
    noise_vars <- setdiff(1:p, true_active)
    info <- dat$variable_info
    
    # ----------------------------
    # Scores raw
    # ----------------------------
    
    sc_raw <- compute_scores(X, Y, M, A)
    
    # ----------------------------
    # Stage 1 : top 10 geo13 raw
    # ----------------------------
    
    ord_geo13_raw_all <- order_geo13(sc_raw, 1:p)
    K1 <- ord_geo13_raw_all[1:k_raw_Y]
    
    remaining1 <- setdiff(1:p, K1)
    
    # Rang Y raw dans remaining1
    ord_Y_raw_R1 <- order_geo13(sc_raw, remaining1)
    rank_Y_raw_R1 <- make_rank(ord_Y_raw_R1, p)
    
    # Rang S2 raw dans remaining1
    ord_S2_R1 <- remaining1[order(-sc_raw$S2[remaining1])]
    rank_S2_R1 <- make_rank(ord_S2_R1, p)
    
    # ----------------------------
    # Stage 2 : top 10 S2 raw
    # ----------------------------
    
    K2 <- ord_S2_R1[1:k_raw_M]
    
    remaining2 <- setdiff(remaining1, K2)
    
    # Rang S2 après avoir retiré K2
    ord_S2_R2 <- remaining2[order(-sc_raw$S2[remaining2])]
    rank_S2_R2 <- make_rank(ord_S2_R2, p)
    
    # ----------------------------
    # Stage 3 : résidualiser Y sur A,M,K1,K2
    # puis top 10 geo13 resid
    # ----------------------------
    
    K12 <- unique(c(K1, K2))
    
    sc_res <- compute_scores_resid_Yonly(X, A, M, Y, K12)
    
    ord_Y_res_R2 <- order_geo13(sc_res, remaining2)
    rank_Y_res_R2 <- make_rank(ord_Y_res_R2, p)
    
    K3 <- ord_Y_res_R2[1:k_res_Y]
    
    Kfinal <- unique(c(K1, K2, K3))
    
    # ----------------------------
    # Stockage
    # ----------------------------
    
    tmp <- data.frame()
    
    for (j in true_active) {
      
      type_j <- info$type[info$variable == j]
      level_j <- info$level[info$variable == j]
      beta_j <- info$beta[info$variable == j]
      
      selected_geo13_raw <- j %in% K1
      selected_s2_raw <- j %in% K2
      selected_geo13_resid <- j %in% K3
      selected_final <- j %in% Kfinal
      
      # Rang outcome Y : seulement si la variable n'est pas déjà prise avant resid
      if (selected_geo13_raw || selected_s2_raw) {
        rank_Y_raw <- NA
        rank_Y_resid <- NA
        rank_Y_gain <- NA
        noise_Y_raw <- NA
        noise_Y_resid <- NA
        noise_Y_gain <- NA
      } else {
        rank_Y_raw <- rank_Y_raw_R1[j]
        rank_Y_resid <- rank_Y_res_R2[j]
        rank_Y_gain <- rank_Y_raw - rank_Y_resid
        
        noise_R2 <- intersect(noise_vars, remaining2)
        
        noise_Y_raw <- sum(rank_Y_raw_R1[noise_R2] < rank_Y_raw, na.rm = TRUE)
        noise_Y_resid <- sum(rank_Y_res_R2[noise_R2] < rank_Y_resid, na.rm = TRUE)
        noise_Y_gain <- noise_Y_raw - noise_Y_resid
      }
      
      # Rang S2 : score naturel pour AM/pM
      if (selected_geo13_raw) {
        rank_S2_before <- NA
        noise_S2_before <- NA
      } else {
        rank_S2_before <- rank_S2_R1[j]
        
        noise_R1 <- intersect(noise_vars, remaining1)
        
        noise_S2_before <- sum(rank_S2_R1[noise_R1] < rank_S2_before,
                               na.rm = TRUE)
      }
      
      tmp <- rbind(
        tmp,
        data.frame(
          run = b,
          variable = j,
          type = type_j,
          level = level_j,
          beta = beta_j,
          is_confounder = j %in% true_conf,
          
          selected_geo13_raw = selected_geo13_raw,
          selected_s2_raw = selected_s2_raw,
          selected_geo13_resid = selected_geo13_resid,
          selected_final = selected_final,
          
          rank_Y_raw = rank_Y_raw,
          rank_Y_resid = rank_Y_resid,
          rank_Y_gain = rank_Y_gain,
          noise_Y_raw = noise_Y_raw,
          noise_Y_resid = noise_Y_resid,
          noise_Y_gain = noise_Y_gain,
          
          rank_S2_before = rank_S2_before,
          noise_S2_before = noise_S2_before
        )
      )
    }
    
    tmp
    
  }, error = function(e) {
    cat("Erreur run", b, ":", e$message, "\n")
    NULL
  })
  
  if (!is.null(out)) {
    res <- rbind(res, out)
  }
}

# ============================================================
# Résultat brut run par run
# ============================================================

print(head(res))