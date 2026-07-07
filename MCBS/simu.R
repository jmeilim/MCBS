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

k_raw_Y <- 10
k_raw_M <- 10
k_res_Y <- 10

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
  
  dat <- DGP_four_signals_corr(n, p, rho_x = 0, noise_sd = 0.3)
  
  X <- dat$X
  A <- dat$A
  M <- dat$M
  Y <- dat$Y
  
  true_active <- dat$true_signals
  info <- dat$variable_info
  
  # garder seulement les signaux liés à Y
  true_Y <- true_active[info$type[match(true_active, info$variable)] %in% c("AY", "MY", "pY")]
  
  # ----------------------------
  # Scores brutes 
  # ----------------------------
  
  sc_raw <- compute_scores(X, Y, M, A)
  
  # ----------------------------
  # Stage 1 : geo13 raw
  # ----------------------------
  
  ord_geo13_raw <- order_geo13(sc_raw, 1:p)
  rank_geo13_raw <- make_rank(ord_geo13_raw, p)
  
  K1 <- ord_geo13_raw[1:k_raw_Y]
  
  remaining1 <- setdiff(1:p, K1)
  
  # --------------------------
  # Stage 2 : S2 raw
  # -------------------------
  
  ord_S2 <- remaining1[order(-sc_raw$S2[remaining1])]
  K2 <- ord_S2[1:k_raw_M]
  
  remaining2 <- setdiff(remaining1, K2)
  
  # --------------------------
  # Stage 3 : geo13 résiduel
  # ----------------------------
  
  K12 <- unique(c(K1, K2))
  
  sc_res <- compute_scores_resid_Yonly(X, A, M, Y, K12)
  
  ord_geo13_resid <- order_geo13(sc_res, remaining2)
  rank_geo13_resid <- make_rank(ord_geo13_resid, p)
  
  # ----------------------------
  # Stockage des rangs 
  # --------------------------
  
  tmp <- data.frame()
  
  for (j in true_Y) {
    
    type_j <- info$type[info$variable == j]
    level_j <- info$level[info$variable == j]
    beta_j <- info$beta[info$variable == j]
    
    selected_stage1 <- j %in% K1
    selected_stage2 <- j %in% K2
    
    rank_raw <- rank_geo13_raw[j]
    
    if (selected_stage1 || selected_stage2) {
      rank_resid <- NA
    } else {
      rank_resid <- rank_geo13_resid[j]
    }
    
    tmp <- rbind(
      tmp,
      data.frame(
        run = b,
        variable = j,
        type = type_j,
        level = level_j,
        beta = beta_j,
        rank_raw = rank_raw,
        rank_resid = rank_resid
      )
    )
  }
  
  tmp
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

summary_ranks <- res %>%
  group_by(variable, type, level, beta) %>%
  summarise(
    mean_rank_raw = mean_na(rank_raw),
    mean_rank_resid = mean_na(rank_resid),
    .groups = "drop"
  ) %>%
  arrange(type, beta)

print(summary_ranks, width = Inf)

summary_ranks <- res %>%
  group_by(variable, type, level, beta) %>%
  summarise(
    mean_rank_raw = median_na(rank_raw),
    mean_rank_resid = median_na(rank_resid),
    .groups = "drop"
  ) %>%
  arrange(type, beta)

print(summary_ranks, width = Inf)