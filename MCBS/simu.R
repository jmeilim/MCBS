library(parallel)
library(dplyr)


resid_lm <- function(y, Z) {
  Z <- as.data.frame(Z)
  fit <- tryCatch(lm(y ~ ., data = data.frame(y = y, Z)),
                  error = function(e) NULL)
  if (is.null(fit)) return(y - mean(y))
  as.numeric(resid(fit))
}

compute_S1_raw <- function(X, Y, A) {
  p <- ncol(X)
  S1 <- numeric(p)
  for (j in 1:p) S1[j] <- Causal.cor(X[, j], Y, A)
  S1
}

compute_S1_res <- function(X, A, Y, K) {
  p <- ncol(X)
  Yres <- resid_lm(Y, cbind(A = A, X[, K, drop = FALSE]))
  S1 <- numeric(p)
  for (j in 1:p) S1[j] <- Causal.cor(X[, j], Yres, A)
  S1
}


order_S1 <- function(S1, candidates) {
  candidates[order(-S1[candidates])]
}

make_rank <- function(ord, p) {
  rk <- rep(NA, p); rk[ord] <- seq_along(ord); rk
}


set.seed(42)

B  <- 100
n  <- 200       
p  <- 1000    
k1 <- 2
k2 <- 15

one_run <- function(b) {
  
  dat <- DGP_cbs_style(n = 200, p = 1000,
                       bA = 2, beta_dom = 5, beta_weak = 1, n_weak = 5,
                       alpha_conf = 2, alpha_instr = 3)
  X <- dat$X; A <- dat$A; Y <- dat$Y
  info <- dat$variable_info
  true_signals <- dat$true_signals
  

  S1_raw <- compute_S1_raw(X, Y, A)
  ord1 <- order_S1(S1_raw, 1:p)
  rk1  <- make_rank(ord1, p)
  K1   <- ord1[1:k1]
  

  remaining <- setdiff(1:p, K1)
  S1_res <- compute_S1_res(X, A, Y, K1)
  ord2   <- order_S1(S1_res, remaining)
  rk2    <- make_rank(ord2, p)
  K2     <- ord2[1:k2]
  
  K_final <- unique(c(K1, K2))
  
  tmp <- data.frame()
  for (j in true_signals) {
    type_j  <- info$type [info$variable == j]
    level_j <- info$level[info$variable == j]
    beta_j  <- info$beta [info$variable == j]
    
    sel1 <- j %in% K1
    sel2 <- j %in% K2
    tmp <- rbind(tmp, data.frame(
      run = b, variable = j,
      type = type_j, level = level_j, beta = beta_j,
      selected_stage1 = sel1,
      selected_stage2 = sel2,
      selected_final  = j %in% K_final,
      rank_stage1 = rk1[j],
      rank_stage2 = ifelse(sel1, NA, rk2[j])
    ))
  }
  tmp
}


n_cores <- 8
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 42)

clusterEvalQ(cl, { .libPaths("~/Rlibs"); library(Ball); library(MASS) })

clusterExport(cl, varlist = c(
  "B", "n", "p", "k1", "k2",
  "DGP_cbs_style", "expit",
  "compute_S1_raw", "compute_S1_res",
  "order_S1", "make_rank", "resid_lm",
  "one_run",
  "Causal.cor"
), envir = environment())

res_list <- parLapply(cl, 1:B, function(b) {
  tryCatch(one_run(b),
           error = function(e) { cat("Erreur run", b, ":", conditionMessage(e), "\n"); NULL })
})
stopCluster(cl)

res <- do.call(rbind, res_list[!sapply(res_list, is.null)])


median_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

summary_median <- res %>%
  group_by(variable, type, level, beta) %>%
  summarise(
    median_rank_stage1 = median_na(rank_stage1),
    median_rank_stage2 = median_na(rank_stage2),
    .groups = "drop"
  ) %>%
  arrange(type, desc(beta))

cat("n =", n, ", p =", p, ", k1 =", k1, ", k2 =", k2, "\n")
print(summary_median, width = Inf)