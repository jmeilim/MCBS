library(pbmcapply)
library(parallel)

N_CORES <- 8

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
  i0 <- which(A == 0); i1 <- which(A == 1)
  yp[i0] <- y[i0][sample(length(i0))]
  yp[i1] <- y[i1][sample(length(i1))]
  yp
}


noise_threshold_par <- function(X, y_cur, A, candidates,
                                B = 200, gamma = 0.10,
                                n_cores = N_CORES) {
  Xc <- X[, candidates, drop = FALSE]   # sous-matrice figée (évite de réindexer dans chaque worker)
  maxima <- pbmclapply(1:B, function(b) {
    yperm <- permute_within_A(y_cur, A)
    sc <- vapply(seq_along(candidates),
                 function(j) Causal.cor(Xc[, j], yperm, A),
                 numeric(1))
    max(sc)
  }, mc.cores = n_cores)
  maxima <- unlist(maxima)
  as.numeric(quantile(maxima, probs = 1 - gamma))
}


scores_par <- function(X, y_cur, A, remaining, n_cores = N_CORES) {
  unlist(pbmclapply(remaining,
                    function(j) Causal.cor(X[, j], y_cur, A),
                    mc.cores = n_cores))
}

peel_until_noise <- function(X, A, Y,
                             gamma       = 0.10,
                             B           = 50,
                             max_iter    = 50,
                             recalibrate = FALSE,
                             n_cores     = N_CORES,
                             verbose     = TRUE) {
  p <- ncol(X)
  peeled    <- integer(0)
  remaining <- 1:p
  y_cur     <- as.numeric(Y)
  
  log_df <- data.frame(iter=integer(0), n_peeled=integer(0),
                       vars_peeled=character(0), score_max=numeric(0),
                       threshold=numeric(0), decision=character(0),
                       stringsAsFactors = FALSE)
  
  # calibration unique au départ (sur tout remaining)
  sub <- sample(remaining, ceiling(1 * length(remaining)))
  thr_fixed <- noise_threshold(X, y_cur, A, sub, B = B, gamma = gamma)
  
  for (it in 1:max_iter) {
    scores <- scores_par(X, y_cur, A, remaining, n_cores = n_cores)
    
    thr <- if (recalibrate) {
      noise_threshold_par(X, y_cur, A, remaining, B = B, gamma = gamma,
                          n_cores = n_cores)
    } else thr_fixed
    
    above_local <- which(scores > thr)
    smax <- max(scores)
    
    if (length(above_local) == 0) {
      log_df <- rbind(log_df, data.frame(iter=it, n_peeled=0L, vars_peeled="",
                                         score_max=smax, threshold=thr, decision="STOP",
                                         stringsAsFactors = FALSE))
      if (verbose) cat(sprintf("iter %2d | paquet vide | score_max=%.3e | seuil=%.3e | STOP\n",
                               it, smax, thr))
      break
    }
    
    vars_this <- remaining[above_local]
    log_df <- rbind(log_df, data.frame(iter=it, n_peeled=length(vars_this),
                                       vars_peeled=paste(vars_this, collapse=","),
                                       score_max=smax, threshold=thr, decision="PEEL",
                                       stringsAsFactors = FALSE))
    if (verbose) cat(sprintf("iter %2d | pele %2d var {%s} | score_max=%.3e | seuil=%.3e | PEEL\n",
                             it, length(vars_this),
                             paste(vars_this, collapse=","), smax, thr))
    
    peeled    <- c(peeled, vars_this)
    remaining <- setdiff(remaining, vars_this)
    y_cur     <- resid_on(as.numeric(Y), A, X, peeled)
    if (length(remaining) == 0) break
  }
  
  list(selected = peeled, k_hat = length(peeled), log = log_df)
}

peel_sequential <- function(X, A, Y,
                            gamma = 0.10, B = 100, max_iter = 50,
                            n_cores = 8, verbose = TRUE) {
  p <- ncol(X)
  peeled <- integer(0); remaining <- 1:p
  y_cur <- as.numeric(Y)
  

  thr <- noise_threshold_par(X, y_cur, A, remaining, B = B, gamma = gamma, n_cores = n_cores)
  if (verbose) cat("Seuil fixe:", format(thr, digits=3), "\n")
  
  for (it in 1:max_iter) {
    scores <- scores_par(X, y_cur, A, remaining, n_cores = n_cores)
    jmax   <- which.max(scores)         
    smax   <- scores[jmax]
    
    if (smax <= thr) {                   
      if (verbose) cat(sprintf("iter %d | max=%.3e <= seuil | STOP\n", it, smax))
      break
    }
    
    var_peeled <- remaining[jmax]
    if (verbose) cat(sprintf("iter %d | pele var %d | score=%.3e > seuil=%.3e\n",
                             it, var_peeled, smax, thr))
    
    peeled    <- c(peeled, var_peeled)
    remaining <- remaining[-jmax]
    y_cur     <- resid_on(as.numeric(Y), A, X, peeled) 
    if (length(remaining) == 0) break
  }
  
  list(selected = peeled, k_hat = length(peeled))
}