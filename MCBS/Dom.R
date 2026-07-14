resid_on <- function(y, A, X, K) {
  if (length(K) == 0) return(y - mean(y))
  Z  <- cbind(A = A, X[, K, drop = FALSE])
  df <- data.frame(y = y, Z)
  fit <- tryCatch(lm(y ~ ., data = df), error = function(e) NULL)
  if (is.null(fit)) return(y - mean(y))
  as.numeric(resid(fit))
}



permute_within_A <- function(y, A) {
  yp <- y
  i0 <- which(A == 0)
  i1 <- which(A == 1)
  yp[i0] <- y[i0][sample(length(i0))]
  yp[i1] <- y[i1][sample(length(i1))]
  yp
}


noise_threshold <- function(X, y_cur, A, candidates,
                            B = 200, gamma = 0.10, verbose = FALSE) {
  maxima <- numeric(B)
  for (b in 1:B) {
    yperm <- permute_within_A(y_cur, A)
    sc <- vapply(candidates,
                 function(j) Causal.cor(X[, j], yperm, A),
                 numeric(1))
    maxima[b] <- max(sc)
  }
  as.numeric(quantile(maxima, probs = 1 - gamma))
}



peel_until_noise <- function(X, A, Y,
                             gamma      = 0.10,   
                             B          = 50,    
                             max_iter   = 50,     
                             recalibrate= FALSE,  
                             verbose    = TRUE ) {
  p <- ncol(X)
  peeled    <- integer(0)     
  remaining <- 1:p
  y_cur     <- as.numeric(Y)
  
  log_df <- data.frame(
    iter        = integer(0),
    n_peeled    = integer(0),   
    vars_peeled = character(0),
    score_max   = numeric(0),
    threshold   = numeric(0),
    decision    = character(0),
    stringsAsFactors = FALSE
  )
  
  sub <- sample(remaining, ceiling(0.50 * length(remaining)))
  thr_fixed <- noise_threshold(X, y_cur, A, sub, B = B, gamma = gamma)
  
  for (it in 1:max_iter) {
    
    scores <- vapply(remaining,
                     function(j) Causal.cor(X[, j], y_cur, A),
                     numeric(1))
    
    thr <- if (recalibrate) {
      noise_threshold(X, y_cur, A, remaining, B = B, gamma = gamma)
    } else {
      thr_fixed
    }
    
    
    above_local <- which(scores > thr)          
    smax        <- max(scores)                 
    
    if (length(above_local) == 0) 
      {
      log_df <- rbind(log_df, data.frame(
        iter = it, n_peeled = 0L, vars_peeled = "",
        score_max = smax, threshold = thr,
        decision = "STOP", stringsAsFactors = FALSE
      ))
      if (verbose) {
        cat(sprintf("iter %2d | paquet vide | score_max=%.3e | seuil=%.3e | STOP\n",
                    it, smax, thr))
      }
      break
    }
    
    vars_this <- remaining[above_local]         
    
    log_df <- rbind(log_df, data.frame(
      iter = it,
      n_peeled = length(vars_this),
      vars_peeled = paste(vars_this, collapse = ","),
      score_max = smax, threshold = thr,
      decision = "PEEL", stringsAsFactors = FALSE
    ))
    
    if (verbose) {
      cat(sprintf("iter %2d | pele %2d var {%s} | score_max=%.3e | seuil=%.3e | PEEL\n",
                  it, length(vars_this),
                  paste(vars_this, collapse = ","), smax, thr))
    }
    
    
    peeled    <- c(peeled, vars_this)
    remaining <- setdiff(remaining, vars_this)
    y_cur     <- resid_on(as.numeric(Y), A, X, peeled)
    
    if (length(remaining) == 0) break           
  }
  
  list(
    selected = peeled,          
    k_hat    = length(peeled),  
    log      = log_df
  )
}


library(Ball)
set.seed(1)
dat <- DGP_cbs_style(n = 200, p = 1000,
                      bA = 2, beta_dom = 7, beta_weak = 2, n_weak = 5,
                      alpha_conf = 2, alpha_instr = 3)

tm <- system.time({
  out <- peel_until_noise(dat$X, dat$A, dat$Y,
                          gamma = 0.10, B = 100, verbose = TRUE)
})

cat("\nk_hat =", out$k_hat, "\n")
cat("Variables pelees :", out$selected, "\n")
cat("Vrais confondeurs:", dat$true_confounders, "\n")
print(out$log)

cat(sprintf("\nTemps total peeling : %.1f s (elapsed)\n", tm["elapsed"]))
