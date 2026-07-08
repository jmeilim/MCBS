# ============================================================
# 10 geo13 raw + 10 S2 raw + 10 geo13 residualized Y
# ============================================================


# ============================================================
# 1. Résidu de Y
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
# 2. Scores après résidualisation de Y 
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
    S2[j] <- Causal.cor(X[, j], M, A)  # M brut
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
# Paramètres Monte Carlo
# ============================================================

set.seed(42)

B <- 100
n <- 100
p <- 1000
lambda <- 4

k1 <- 5
k2 <- 15

res <- data.frame()

# ============================================================
# Moyenne avec NA 
# ============================================================

mean_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(mean(x, na.rm = TRUE))
  }
}

median_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(median(x, na.rm = TRUE))
  }
}

# ============================================================
# Une seule run Monte Carlo
# ============================================================

one_run <- function(b) {
  
  dat <- DGP_four_signals_corr(n, p, rho_x = 0, noise_sd = 0.3 )
  
  X <- dat$X
  A <- dat$A
  M <- dat$M
  Y <- dat$Y
  
  info <- dat$variable_info
  true_signals <- dat$true_signals
  
  # ============================================================
  # Stage 1 : geo13 raw
  # ============================================================
  
  sc_raw <- compute_scores(X, Y, M, A)
  
  ord_stage1 <- order_geo13(sc_raw, 1:p)
  rank_stage1 <- make_rank(ord_stage1, p)
  
  K1 <- ord_stage1[1:k1]
  
  # ============================================================
  # Stage 2 : geo13 après résidualisation sur K1
  # ============================================================
  
  remaining1 <- setdiff(1:p, K1)
  
  sc_res <- compute_scores_resid_Yonly(X, A, M, Y, K1)
  
  ord_stage2 <- order_geo13(sc_res, remaining1)
  rank_stage2 <- make_rank(ord_stage2, p)
  
  K2 <- ord_stage2[1:k2]
  
  # ============================================================
  # Ensemble final
  # ============================================================
  
  K_final <- unique(c(K1, K2))
  
  # ============================================================
  # Résultats pour les vraies variables
  # ============================================================
  
  tmp <- data.frame()
  
  for (j in true_signals) {
    
    type_j <- info$type[info$variable == j]
    level_j <- info$level[info$variable == j]
    beta_j <- info$beta[info$variable == j]
    
    selected_stage1 <- j %in% K1
    selected_stage2 <- j %in% K2
    selected_final <- j %in% K_final
    
    tmp <- rbind(
      tmp,
      data.frame(
        run = b,
        variable = j,
        type = type_j,
        level = level_j,
        beta = beta_j,
        
        selected_stage1 = selected_stage1,
        selected_stage2 = selected_stage2,
        selected_final = selected_final,
        
        rank_stage1 = rank_stage1[j],
        rank_stage2 = ifelse(selected_stage1, NA, rank_stage2[j])
      )
    )
  }
  
  return(tmp)
}

# ============================================================
# Boucle Monte Carlo
# ============================================================

for (b in 1:B) {
  
  cat("Run", b, "/", B, "\n")
  
  out <- tryCatch(
    one_run(b),
    error = function(e) {
      cat("Erreur run", b, ":", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(out)) {
    res <- rbind(res, out)
  }
}

# ============================================================
# Résumé : rangs moyens 
# ============================================================

library(dplyr)

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

summary_mean <- res %>%
  group_by(variable, type, level, beta) %>%
  summarise(
    mean_rank_stage1 = mean_na(rank_stage1),
    mean_rank_stage2 = mean_na(rank_stage2),
    
    p_stage1 = mean(selected_stage1),
    p_stage2 = mean(selected_stage2),
    p_final = mean(selected_final),
    
    .groups = "drop"
  ) %>%
  arrange(type, beta)

print(summary_mean, width = Inf)
median_na <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

summary_median <- res %>%
  group_by(variable, type, level, beta) %>%
  summarise(
    median_rank_stage1 = median_na(rank_stage1),
    median_rank_stage2 = median_na(rank_stage2),
    
    p_stage1 = mean(selected_stage1),
    p_stage2 = mean(selected_stage2),
    p_final = mean(selected_final),
    
    .groups = "drop"
  ) %>%
  arrange(type, beta)

print(summary_median, width = Inf)
