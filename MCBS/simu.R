source("MCBS/MCBS.R")
source("MCBS/DGP.R")  

# ============================================================
# COMPARAISON RAPIDE CIBLÉE
# raw30 vs hybrid 20-10 vs hybrid 15-10-5
# n = 300, p = 1000, lambda = 4, B = 100, k = 30
# ============================================================

# ------------------------------------------------------------
# 0. Utilitaires
# ------------------------------------------------------------



safe_lm_resid <- function(y, Z) {
  Z <- as.data.frame(Z)
  df <- data.frame(y = y, Z)
  
  fit <- tryCatch(
    lm(y ~ ., data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(y - mean(y))
  }
  
  resid(fit)
}

safe_glm_resid <- function(y, Z) {
  Z <- as.data.frame(Z)
  df <- data.frame(y = y, Z)
  
  fit <- tryCatch(
    suppressWarnings(glm(y ~ ., data = df, family = binomial())),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(y - mean(y))
  }
  
  phat <- fitted(fit)
  phat <- pmin(pmax(phat, 1e-4), 1 - 1e-4)
  
  y - phat
}

# ------------------------------------------------------------
# 1. Scores résidualisés depuis K0
# ------------------------------------------------------------

compute_scores_resid_from_K0 <- function(X, Y, M, A, K0) {
  p <- ncol(X)
  
  if (length(K0) == 0) {
    Y_res <- Y
    M_res <- M
  } else {
    ZY <- cbind(A = A, M = M, X[, K0, drop = FALSE])
    ZM <- cbind(A = A, X[, K0, drop = FALSE])
    
    Y_res <- safe_lm_resid(Y, ZY)
    M_res <- safe_glm_resid(M, ZM)
  }
  
  s1 <- s2 <- s3 <- numeric(p)
  
  for (j in 1:p) {
    s1[j] <- Causal.cor(X[, j], Y_res, A)
    s2[j] <- Causal.cor(X[, j], M_res, A)
    s3[j] <- MCBS.cor(as.matrix(X[, j]), as.matrix(Y_res), A, M)
  }
  
  list(S1 = s1, S2 = s2, S3 = s3)
}

# ------------------------------------------------------------
# 2. Sélection sur sous-ensemble de variables
# ------------------------------------------------------------

select_from_scores_subset <- function(sc, remaining, k, method) {
  if (k <= 0) return(integer(0))
  
  sc_sub <- list(
    S1 = sc$S1[remaining],
    S2 = sc$S2[remaining],
    S3 = sc$S3[remaining]
  )
  
  if (method == "geo") {
    local <- screening_geo(sc_sub, k = k)
  } else if (method == "deux_etapes") {
    local <- screening_deux_etapes(sc_sub, p = length(remaining), k = k)
  } else if (method == "maxZ") {
    local <- screening_final_maxZ(sc_sub, k = k)
  } else {
    stop("Unknown method: ", method)
  }
  
  remaining[local]
}

select_S2_resid_subset <- function(sc_resid, remaining, k) {
  if (k <= 0) return(integer(0))
  
  s2 <- sc_resid$S2[remaining]
  remaining[order(-s2)[1:k]]
}

# ------------------------------------------------------------
# 3. Évaluation globale
# ------------------------------------------------------------

eval_selection_global <- function(selected, dgp) {
  info <- dgp$variable_info
  
  conf <- dgp$true_confounders
  active <- dgp$true_signals
  
  AY <- info$variable[info$type == "AY"]
  AM <- info$variable[info$type == "AM"]
  MY <- info$variable[info$type == "MY"]
  pY <- info$variable[info$type == "pY"]
  pM <- info$variable[info$type == "pM"]
  
  weak_active <- info$variable[info$level == "weak"]
  medium_active <- info$variable[info$level == "medium"]
  strong_active <- info$variable[info$level == "strong"]
  
  weak_conf <- info$variable[
    info$level == "weak" & info$type %in% c("AY", "AM", "MY")
  ]
  
  data.frame(
    n_active = sum(active %in% selected),
    n_confounders = sum(conf %in% selected),
    
    prop_active = sum(active %in% selected) / length(active),
    prop_confounders = sum(conf %in% selected) / length(conf),
    
    prop_AY = sum(AY %in% selected) / length(AY),
    prop_AM = sum(AM %in% selected) / length(AM),
    prop_MY = sum(MY %in% selected) / length(MY),
    prop_pY = sum(pY %in% selected) / length(pY),
    prop_pM = sum(pM %in% selected) / length(pM),
    
    prop_weak_active = sum(weak_active %in% selected) / length(weak_active),
    prop_medium_active = sum(medium_active %in% selected) / length(medium_active),
    prop_strong_active = sum(strong_active %in% selected) / length(strong_active),
    
    prop_weak_conf = sum(weak_conf %in% selected) / length(weak_conf),
    
    all_active = as.numeric(all(active %in% selected)),
    all_confounders = as.numeric(all(conf %in% selected))
  )
}

eval_selection_by_var <- function(selected, dgp) {
  info <- dgp$variable_info
  
  data.frame(
    variable = info$variable,
    type = info$type,
    level = info$level,
    beta = info$beta,
    selected = as.numeric(info$variable %in% selected)
  )
}

# ------------------------------------------------------------
# 4. Simulation rapide ciblée
# ------------------------------------------------------------

compare_fast_targeted <- function(
    B = 100,
    n = 600,
    p = 2000,
    lambda = 4,
    seed = 123,
    verbose = TRUE
) {
  set.seed(seed)
  
  global_results <- list()
  var_results <- list()
  
  for (b in 1:B) {
    if (verbose) {
      cat("Run", b, "/", B, "\n")
    }
    
    dgp <- DGP_rich(n = n, p = p, lambda = lambda)
    
    # Scores bruts calculés une seule fois
    sc_raw <- compute_scores(dgp$X, dgp$Y, dgp$M, dgp$A)
    
    selections <- list()
    
    # --------------------------------------------------------
    # A. Baselines raw k = 30
    # --------------------------------------------------------
    
    selections[["raw30_deux_etapes"]] <- screening_deux_etapes(sc_raw, p = p, k = 30)
    selections[["raw30_geo"]] <- screening_geo(sc_raw, k = 30)
    selections[["raw30_maxZ"]] <- screening_final_maxZ(sc_raw, k = 30)
    
    # --------------------------------------------------------
    # B. Hybrid 20 + 10
    # raw = deux_etapes
    # resid = maxZ / deux_etapes
    # --------------------------------------------------------
    
    K20 <- screening_deux_etapes(sc_raw, p = p, k = 20)
    
    sc_resid20 <- compute_scores_resid_from_K0(
      X = dgp$X,
      Y = dgp$Y,
      M = dgp$M,
      A = dgp$A,
      K0 = K20
    )
    
    remaining20 <- setdiff(1:p, K20)
    
    Kadd20_maxZ <- select_from_scores_subset(
      sc = sc_resid20,
      remaining = remaining20,
      k = 10,
      method = "maxZ"
    )
    
    selections[["hyb20_10_deux_maxZ"]] <- unique(c(K20, Kadd20_maxZ))[1:30]
    
    Kadd20_deux <- select_from_scores_subset(
      sc = sc_resid20,
      remaining = remaining20,
      k = 10,
      method = "deux_etapes"
    )
    
    selections[["hyb20_10_deux_deux"]] <- unique(c(K20, Kadd20_deux))[1:30]
    
    # --------------------------------------------------------
    # C. Hybrid 15 + 10 + 5
    # raw = deux_etapes
    # yres = maxZ / deux_etapes
    # mres = S2_resid
    # --------------------------------------------------------
    
    K15 <- screening_deux_etapes(sc_raw, p = p, k = 15)
    
    sc_resid15 <- compute_scores_resid_from_K0(
      X = dgp$X,
      Y = dgp$Y,
      M = dgp$M,
      A = dgp$A,
      K0 = K15
    )
    
    remaining15 <- setdiff(1:p, K15)
    
    # 15 + 10 maxZ + 5 S2
    K_yres_maxZ <- select_from_scores_subset(
      sc = sc_resid15,
      remaining = remaining15,
      k = 10,
      method = "maxZ"
    )
    
    remaining15_maxZ <- setdiff(remaining15, K_yres_maxZ)
    
    K_mres_maxZ <- select_S2_resid_subset(
      sc_resid = sc_resid15,
      remaining = remaining15_maxZ,
      k = 5
    )
    
    selections[["hyb15_10_5_deux_maxZ_S2"]] <- unique(
      c(K15, K_yres_maxZ, K_mres_maxZ)
    )[1:30]
    
    # 15 + 10 deux_etapes + 5 S2
    K_yres_deux <- select_from_scores_subset(
      sc = sc_resid15,
      remaining = remaining15,
      k = 10,
      method = "deux_etapes"
    )
    
    remaining15_deux <- setdiff(remaining15, K_yres_deux)
    
    K_mres_deux <- select_S2_resid_subset(
      sc_resid = sc_resid15,
      remaining = remaining15_deux,
      k = 5
    )
    
    selections[["hyb15_10_5_deux_deux_S2"]] <- unique(
      c(K15, K_yres_deux, K_mres_deux)
    )[1:30]
    
    # --------------------------------------------------------
    # Évaluation
    # --------------------------------------------------------
    
    for (m in names(selections)) {
      global_results[[length(global_results) + 1]] <- cbind(
        run = b,
        method = m,
        lambda = lambda,
        k = 30,
        eval_selection_global(selections[[m]], dgp)
      )
      
      var_results[[length(var_results) + 1]] <- cbind(
        run = b,
        method = m,
        lambda = lambda,
        k = 30,
        eval_selection_by_var(selections[[m]], dgp)
      )
    }
  }
  
  list(
    global = do.call(rbind, global_results),
    by_var = do.call(rbind, var_results)
  )
}

# ------------------------------------------------------------
# 5. Summaries
# ------------------------------------------------------------

make_fast_summary <- function(res) {
  main_summary <- aggregate(
    cbind(
      prop_active,
      prop_confounders,
      prop_AY,
      prop_AM,
      prop_MY,
      prop_pY,
      prop_pM,
      prop_weak_active,
      prop_medium_active,
      prop_strong_active,
      prop_weak_conf,
      n_active,
      n_confounders
    ) ~ method,
    data = res$global,
    FUN = mean
  )
  
  main_summary <- main_summary[order(-main_summary$prop_active), ]
  
  type_summary <- aggregate(
    selected ~ method + type,
    data = res$by_var,
    FUN = mean
  )
  
  type_summary <- type_summary[order(type_summary$type, -type_summary$selected), ]
  
  level_summary <- aggregate(
    selected ~ method + level,
    data = res$by_var,
    FUN = mean
  )
  
  level_summary <- level_summary[order(level_summary$level, -level_summary$selected), ]
  
  type_level_summary <- aggregate(
    selected ~ method + type + level,
    data = res$by_var,
    FUN = mean
  )
  
  weak_conf_summary <- subset(
    type_level_summary,
    type %in% c("AY", "AM", "MY") & level == "weak"
  )
  
  weak_conf_summary <- weak_conf_summary[
    order(weak_conf_summary$type, -weak_conf_summary$selected),
  ]
  
  list(
    main = main_summary,
    type = type_summary,
    level = level_summary,
    type_level = type_level_summary,
    weak_conf = weak_conf_summary
  )
}

round_numeric_df <- function(df, digits = 3) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
  df
}

print_fast_summary <- function(summ, digits = 3) {
  cat("\n============================================================\n")
  cat("1. Résumé principal\n")
  cat("============================================================\n")
  print(round_numeric_df(summ$main, digits), row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("2. Résumé par type\n")
  cat("============================================================\n")
  print(round_numeric_df(summ$type, digits), row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("3. Résumé par niveau\n")
  cat("============================================================\n")
  print(round_numeric_df(summ$level, digits), row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("4. Focus — confondeurs faibles\n")
  cat("============================================================\n")
  print(round_numeric_df(summ$weak_conf, digits), row.names = FALSE)
}

# ------------------------------------------------------------
# 6. Lancement
# ------------------------------------------------------------

res_fast <- compare_fast_targeted(
  B = 100,
  n = 600,
  p = 2000,
  lambda = 4,
  seed = 123,
  verbose = TRUE
)

summ_fast <- make_fast_summary(res_fast)

print_fast_summary(summ_fast, digits = 3)