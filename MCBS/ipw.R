
ipw_from_selected <- function(X, D, Y, K, trim = 0.02, hajek = TRUE) {
  if (length(K) == 0) return(NA_real_)
  Xs <- as.data.frame(X[, K, drop = FALSE])
  names(Xs) <- paste0("v", seq_along(K))
  
  psfit <- suppressWarnings(
    glm(D ~ ., data = data.frame(D = D, Xs), family = binomial))
  ps <- as.numeric(predict(psfit, type = "response"))
  ps <- pmin(pmax(ps, trim), 1 - trim)
  
  w <- create_weights(ps, D)
  
  if (hajek) {
    m1 <- sum(w[D == 1] * Y[D == 1]) / sum(w[D == 1])
    m0 <- sum(w[D == 0] * Y[D == 0]) / sum(w[D == 0])
  } else {
    m1 <- sum(w[D == 1] * Y[D == 1]) / length(Y)
    m0 <- sum(w[D == 0] * Y[D == 0]) / length(Y)
  }
  m1 - m0
}

dr_from_selected <- function(X, D, Y, K, trim = 0.02) {
  if (length(K) == 0) return(NA_real_)
  Xs <- as.data.frame(X[, K, drop = FALSE])
  names(Xs) <- paste0("v", seq_along(K))
  
  psfit <- suppressWarnings(
    glm(D ~ ., data = data.frame(D = D, Xs), family = binomial))
  ps <- as.numeric(predict(psfit, type = "response"))
  ps <- pmin(pmax(ps, trim), 1 - trim)
  
  fit_or <- function(arm_idx) {
    dtr <- data.frame(Y = Y[arm_idx], Xs[arm_idx, , drop = FALSE])
    m   <- lm(Y ~ ., data = dtr)
    as.numeric(predict(m, newdata = Xs))    
  }
  b1 <- fit_or(which(D == 1))
  b0 <- fit_or(which(D == 0))
  
  mean(DR(ps, b1, b0, D, Y))
}

run_dr <- function(seed, n = 300, p = 1000, ac = 1.0, bd = 10, bc = 2) {
  set.seed(seed)
  X <- matrix(runif(n*p, -1, 1), ncol = p)
  ps_true <- 1/(1 + exp(-rowSums(X[, 1:2])*ac))
  D <- rbinom(n, 1, ps_true)
  Y <- D*2 + rowSums(X[, 1:2])*bc + X[, 4]*bd + rnorm(n)
  
  target <- c(1, 2, 4)
  
  # screening CBS : top-30
  rho   <- vapply(1:p, function(j) Causal.cor(X[, j], Y, D), numeric(1))
  K_cbs <- order(rho, decreasing = TRUE)[1:min(30, p)]
  
  # peeling
  out    <- peel_until_noise(X, D, Y, gamma = 0.10, B = 100, verbose = FALSE)
  K_peel <- out$selected
  
  c(dr_oracle       = dr_from_selected(X, D, Y, c(1, 2)),
    dr_oracle_plus  = dr_from_selected(X, D, Y, target),
    dr_cbs          = dr_from_selected(X, D, Y, K_cbs),
    dr_peel         = dr_from_selected(X, D, Y, K_peel),
    ipw_oracle      = ipw_from_selected(X, D, Y, c(1, 2)),
    ipw_oracle_plus = ipw_from_selected(X, D, Y, target),
    ipw_cbs         = ipw_from_selected(X, D, Y, K_cbs),
    ipw_peel        = ipw_from_selected(X, D, Y, K_peel),
    k_hat           = out$k_hat,
    cap_cbs         = as.numeric(all(c(1, 2) %in% K_cbs)),
    cap_peel        = as.numeric(all(c(1, 2) %in% K_peel)),
    exact_target    = as.numeric(setequal(K_peel, target)))
}


library(parallel)
run_mc_dr <- function(seeds, ncores = min(detectCores() - 1, 8)) {
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, { library(Ball); library(glmnet) })
  clusterExport(cl, c("run_dr", "dr_from_selected", "ipw_from_selected",
                      "create_weights", "peel_until_noise",
                      "noise_threshold", "permute_within_A", "resid_on",
                      "Causal.cor", "DR", "expit"), envir = .GlobalEnv)
  t(parSapply(cl, seeds, run_dr))
}


MC2 <- as.data.frame(run_mc_dr(1:100))
truth <- 2
cat(sprintf("beta_dom=10  (%d runs)\n\n", nrow(MC2)))
res2 <- t(sapply(c("dr_oracle","dr_oracle_plus","dr_cbs","dr_peel",
                   "ipw_oracle","ipw_oracle_plus","ipw_cbs","ipw_peel"),
                 function(nm){ v <- MC2[[nm]]
                 c(biais = mean(v)-truth, mc_se = sd(v)/sqrt(length(v)),
                   sd = sd(v), rmse = sqrt(mean((v-truth)^2))) }))
print(round(res2, 4))
cat(sprintf("\ncapture {1,2}  CBS %.0f%%  |  Peeling %.0f%%\n",
            100*mean(MC2$cap_cbs), 100*mean(MC2$cap_peel)))
cat(sprintf("cible exacte {1,2,4} : %.0f%%   | k_hat median %.0f [%.0f,%.0f]\n",
            100*mean(MC2$exact_target),
            median(MC2$k_hat), min(MC2$k_hat), max(MC2$k_hat)))