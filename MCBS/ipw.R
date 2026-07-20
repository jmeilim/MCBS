# ============================================================
# Piste 1 : IPW pur — 100 runs
# L'estimateur ne depend QUE du PS, donc uniquement de la selection.
# ============================================================

ipw_from_selected <- function(X, D, Y, K, trim = 0.02) {
  if (length(K) == 0) return(NA_real_)
  Xs <- as.data.frame(X[, K, drop = FALSE])
  names(Xs) <- paste0("v", seq_along(K))
  fit <- suppressWarnings(
    glm(D ~ ., data = data.frame(D = D, Xs), family = binomial))
  ps <- as.numeric(predict(fit, type = "response"))
  ps <- pmin(pmax(ps, trim), 1 - trim)
  mean(D*Y/ps) - mean((1-D)*Y/(1-ps))
}

run_ipw <- function(seed, n = 200, p = 1000, ac = 1.0, bd = 10, bc = 2) {
  set.seed(seed)
  X <- matrix(runif(n*p, -1, 1), ncol = p)
  ps_true <- 1/(1 + exp(-rowSums(X[, 1:2])*ac))
  D <- rbinom(n, 1, ps_true)
  Y <- D*2 + rowSums(X[, 1:2])*bc + X[, 4]*bd + rnorm(n)
  
  rho   <- vapply(1:p, function(j) Causal.cor(X[, j], Y, D), numeric(1))
  K_cbs <- order(rho, decreasing = TRUE)[1:min(30, p)]
  out   <- peel_until_noise(X, D, Y, gamma = 0.10, B = 100, verbose = FALSE)
  
  c(ipw_pstrue      = mean(D*Y/ps_true) - mean((1-D)*Y/(1-ps_true)),
    ipw_oracle      = ipw_from_selected(X, D, Y, c(1, 2)),
    ipw_oracle_plus = ipw_from_selected(X, D, Y, c(1, 2, 4)),
    ipw_cbs         = ipw_from_selected(X, D, Y, K_cbs),
    ipw_peel        = ipw_from_selected(X, D, Y, out$selected),
    k_hat           = out$k_hat,
    cap_cbs         = as.numeric(all(c(1, 2) %in% K_cbs)),
    cap_peel        = as.numeric(all(c(1, 2) %in% out$selected)),
    has_X4          = as.numeric(4 %in% out$selected),
    exact_target    = as.numeric(setequal(out$selected, c(1, 2, 4))))
}

# --- Monte Carlo (PSOCK, Windows-compatible) ---
library(parallel)
run_mc_ipw <- function(seeds, ncores = min(detectCores() - 1, 8)) {
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, { library(Ball); library(glmnet) })
  clusterExport(cl, c("run_ipw", "ipw_from_selected", "peel_until_noise",
                      "noise_threshold", "permute_within_A", "resid_on",
                      "Causal.cor", "expit"), envir = .GlobalEnv)
  t(parSapply(cl, seeds, run_ipw))
}

MC <- as.data.frame(run_mc_ipw(1:200))

# --- resume ---
truth <- 2
cat(sprintf("IPW  n=200  p=1000  alpha_conf=1.0  beta_dom=10  (%d runs)\n\n", nrow(MC)))
cat(sprintf("%-16s %9s %8s %8s\n", "estimateur", "biais", "sd", "rmse"))
for (nm in c("ipw_pstrue", "ipw_oracle", "ipw_oracle_plus", "ipw_cbs", "ipw_peel")) {
  v <- MC[[nm]]
  cat(sprintf("%-16s %+9.4f %8.4f %8.4f\n",
              nm, mean(v) - truth, sd(v), sqrt(mean((v - truth)^2))))
}
cat(sprintf("\ncapture {1,2}   CBS %3.0f%%  |  Peeling %3.0f%%\n",
            100*mean(MC$cap_cbs), 100*mean(MC$cap_peel)))
cat(sprintf("X4 (precision) dans le pele : %3.0f%%\n", 100*mean(MC$has_X4)))
cat(sprintf("cible exacte {1,2,4}        : %3.0f%%\n", 100*mean(MC$exact_target)))
cat(sprintf("k_hat : median %.0f  [%.0f, %.0f]\n",
            median(MC$k_hat), min(MC$k_hat), max(MC$k_hat)))