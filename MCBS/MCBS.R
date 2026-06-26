# ============================================================
# MCBS.R — CBS Mediation
# ============================================================
library(Ball)
library(glmnet)


MCBS.cor <- function(x, y, A, M, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  idx00 <- which(A==0 & M==0); idx01 <- which(A==0 & M==1)
  idx10 <- which(A==1 & M==0); idx11 <- which(A==1 & M==1)
  p00 <- length(idx00)/length(A); p01 <- length(idx01)/length(A)
  p10 <- length(idx10)/length(A); p11 <- length(idx11)/length(A)
  if (distance == TRUE) {
    MCBS.cor <-
      p00*bcov(x[idx00,idx00], y[idx00,idx00], distance=TRUE)^2 +
      p01*bcov(x[idx01,idx01], y[idx01,idx01], distance=TRUE)^2 +
      p10*bcov(x[idx10,idx10], y[idx10,idx10], distance=TRUE)^2 +
      p11*bcov(x[idx11,idx11], y[idx11,idx11], distance=TRUE)^2
  } else {
    MCBS.cor <-
      p00*bcov(x[idx00,], y[idx00,])^2 +
      p01*bcov(x[idx01,], y[idx01,])^2 +
      p10*bcov(x[idx10,], y[idx10,])^2 +
      p11*bcov(x[idx11,], y[idx11,])^2
  }
  return(MCBS.cor)
}


Causal.cor <- function(x, y, z, distance = FALSE) {
  x      <- as.matrix(x)
  y      <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha  <- length(index0) / length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]; y0 <- y[index0, index0]
    x1 <- x[index1, index1]; y1 <- y[index1, index1]
    Causal.cor <- alpha*(1-alpha)*bcov(x0,y0,distance=TRUE)^2 +
      (1-alpha)*alpha*bcov(x1,y1,distance=TRUE)^2
  } else {
    x0 <- x[index0, ]; y0 <- y[index0, ]
    x1 <- x[index1, ]; y1 <- y[index1, ]
    Causal.cor <- alpha*bcov(x0,y0)^2 + (1-alpha)*bcov(x1,y1)^2
  }
  return(Causal.cor)
}


# ── Screening ───────────────────────────────────────────────

compute_scores <- function(X, Y, M, A) {
  p  <- ncol(X)
  s1 <- s2 <- s3 <- numeric(p)
  for(j in 1:p) {
    s1[j] <- Causal.cor(X[,j], Y, A)
    s2[j] <- Causal.cor(X[,j], M, A)
    s3[j] <- MCBS.cor(as.matrix(X[,j]), as.matrix(Y), A, M)
  }
  list(S1=s1, S2=s2, S3=s3)
}

screening_geo <- function(sc, k) {
  order((rank(-sc$S1)*rank(-sc$S2)*rank(-sc$S3))^(1/3))[1:k]
}

screening_deux_etapes <- function(sc, p, k) {
  k2  <- floor(k/3); k1 <- k - k2
  top <- order(-sc$S2)[1:k2]
  rst <- setdiff(1:p, top)
  geo <- order((rank(-sc$S1[rst])*rank(-sc$S3[rst]))^0.5)[1:k1]
  union(top, rst[geo])
}

# ------------------------------------------------------------
# max-Z
# ------------------------------------------------------------

robust_z <- function(s) {
  med <- median(s, na.rm = TRUE)
  sc  <- mad(s, constant = 1, na.rm = TRUE)
  
  if (!is.finite(sc) || sc < 1e-12) {
    sc <- sd(s, na.rm = TRUE)
  }
  if (!is.finite(sc) || sc < 1e-12) {
    sc <- 1
  }
  
  (s - med) / sc
}

screening_final_maxZ <- function(sc, k) {
  z1 <- robust_z(sc$S1)
  z2 <- robust_z(sc$S2)
  z3 <- robust_z(sc$S3)
  
  score <- pmax(z1, z2, z3)
  order(-score)[1:k]
}


# ------------------------------------------------------------
# Min-rank
# ------------------------------------------------------------

screening_min_rank <- function(sc, k) {
  r1 <- rank(-sc$S1, ties.method = "average")
  r2 <- rank(-sc$S2, ties.method = "average")
  r3 <- rank(-sc$S3, ties.method = "average")
  
  rmin <- pmin(r1, r2, r3)
  order(rmin)[1:k]
}


# ------------------------------------------------------------
# Union blocks
# Par défaut : 30% S1, 40% S2, 30% S3
# ------------------------------------------------------------

screening_union_blocks <- function(
    sc,
    k,
    w = c(S1 = 0.30, S2 = 0.40, S3 = 0.30)
) {
  p <- length(sc$S1)
  
  k1 <- floor(w["S1"] * k)
  k2 <- floor(w["S2"] * k)
  k3 <- k - k1 - k2
  
  K1 <- order(-sc$S1)[1:k1]
  K2 <- order(-sc$S2)[1:k2]
  K3 <- order(-sc$S3)[1:k3]
  
  K <- unique(c(K1, K2, K3))
  
  # Rang minimal pour compléter ou couper
  r1 <- rank(-sc$S1, ties.method = "average")
  r2 <- rank(-sc$S2, ties.method = "average")
  r3 <- rank(-sc$S3, ties.method = "average")
  rmin <- pmin(r1, r2, r3)
  
  # Si doublons => l'union est plus petite que k, on complète
  if (length(K) < k) {
    remaining <- setdiff(1:p, K)
    add <- remaining[order(rmin[remaining])][1:(k - length(K))]
    K <- c(K, add)
  }
  
  # Si jamais l'union dépasse k, on coupe par min-rank
  if (length(K) > k) {
    K <- K[order(rmin[K])][1:k]
  }
  
  K
}

# ── Estimation ──────────────────────────────────────────────

DR <- function(ps, b1, b0, A, Y)
  (A*Y-(A-ps)*b1)/ps - ((1-A)*Y+(A-ps)*b0)/(1-ps)

oal <- function(X, A, w) {
  n  <- nrow(X); vl <- paste0("V",1:ncol(X))
  Xs <- as.data.frame(scale(X)); names(Xs) <- vl; Xs$A <- A
  lv <- seq(log(max(ncol(X),n))^.75/n^.5*.1,
            log(max(ncol(X),n))^.75/n^.5*10, l=10)
  best_ps <- NULL; best_wAMD <- Inf
  for(il in lv) for(ig in 3:20) {
    ps <- pmax(pmin(as.vector(predict(
      glmnet(as.matrix(Xs[,vl]),A,family="binomial",
             alpha=1,penalty.factor=abs(w)^(-ig),lambda=il),
      as.matrix(Xs[,vl]),type="response")),0.95),0.05)
    wt <- ifelse(A==1,1/ps,1/(1-ps))
    Xs$wt <- wt
    wAMD <- sum(abs(w)*sapply(vl,function(v)
      abs(sum(Xs[A==1,v]*wt[A==1])/sum(wt[A==1]) -
            sum(Xs[A==0,v]*wt[A==0])/sum(wt[A==0]))))
    if(wAMD < best_wAMD) { best_wAMD <- wAMD; best_ps <- ps }
  }
  best_ps
}

dr_est <- function(ps, X, A, Y) {
  or <- function(Xtr, Ytr) {
    cv  <- cv.glmnet(Xtr, Ytr, alpha=1, nfold=min(10,nrow(Xtr)))
    idx <- which(coef(cv,"lambda.min")[-1]!=0)
    if(!length(idx)) return(rep(mean(Ytr),nrow(X)))
    predict(lm(Ytr~Xtr[,idx]), newdata=data.frame(X[,idx]))
  }
  dr <- DR(ps, or(X[A==1,],Y[A==1]), or(X[A==0,],Y[A==0]), A, Y)
  c(est=mean(dr), se=sqrt(var(dr)/length(Y)))
}

# ── Pipeline principal ───────────────────────────────────────

CBS_mediation_compare_resid <- function(X, A, M, Y, k = 30, alpha = 0.05) {
  
  p <- ncol(X)
  
  k <- min(k, p)
  
  
  
  nm <- function(s) s / (max(s) + 1e-10)
  
  z  <- qnorm(1 - alpha / 2)
  
  
  
  # Scores bruts
  
  sc <- compute_scores(X, Y, M, A)
  
  
  
  # ----------------------------------------------------------
  
  # 1. Selections brutes
  
  # ----------------------------------------------------------
  
  
  
  sel_2step <- screening_deux_etapes(sc, p, k)
  
  sel_geo   <- screening_geo(sc, k)
  
  sel_maxZ  <- screening_final_maxZ(sc, k)
  
  
  
  # ----------------------------------------------------------
  
  # 2. Selection hybride 20 + 10
  
  # raw = deux_etapes, resid = deux_etapes
  
  # ----------------------------------------------------------
  
  
  
  K20 <- screening_deux_etapes(sc, p = p, k = 20)
  
  
  
  sc_resid20 <- compute_scores_resid_from_K0(
    
    X = X,
    
    Y = Y,
    
    M = M,
    
    A = A,
    
    K0 = K20
    
  )
  
  
  
  remaining20 <- setdiff(1:p, K20)
  
  
  
  Kadd20_deux <- select_from_scores_subset(
    
    sc = sc_resid20,
    
    remaining = remaining20,
    
    k = 10,
    
    method = "deux_etapes"
    
  )
  
  
  
  sel_hyb20_10_deux <- unique(c(K20, Kadd20_deux))[1:k]
  
  
  
  # ----------------------------------------------------------
  
  # 3. Selection hybride 20 + 10
  
  # raw = deux_etapes, resid = maxZ
  
  # ----------------------------------------------------------
  
  
  
  Kadd20_maxZ <- select_from_scores_subset(
    
    sc = sc_resid20,
    
    remaining = remaining20,
    
    k = 10,
    
    method = "maxZ"
    
  )
  
  
  
  sel_hyb20_10_maxZ <- unique(c(K20, Kadd20_maxZ))[1:k]
  
  
  
  # ----------------------------------------------------------
  
  # 4. Selection hybride 15 + 10 + 5
  
  # raw = deux_etapes, resid = deux_etapes, S2_resid
  
  # ----------------------------------------------------------
  
  
  
  K15 <- screening_deux_etapes(sc, p = p, k = 15)
  
  
  
  sc_resid15 <- compute_scores_resid_from_K0(
    
    X = X,
    
    Y = Y,
    
    M = M,
    
    A = A,
    
    K0 = K15
    
  )
  
  
  
  remaining15 <- setdiff(1:p, K15)
  
  
  
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
  
  
  
  sel_hyb15_10_5_deux_S2 <- unique(
    
    c(K15, K_yres_deux, K_mres_deux)
    
  )[1:k]
  
  
  
  # ----------------------------------------------------------
  
  # 5. Selection hybride 15 + 10 + 5
  
  # raw = deux_etapes, resid = maxZ, S2_resid
  
  # ----------------------------------------------------------
  
  
  
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
  
  
  
  sel_hyb15_10_5_maxZ_S2 <- unique(
    
    c(K15, K_yres_maxZ, K_mres_maxZ)
    
  )[1:k]
  
  
  
  # ----------------------------------------------------------
  
  # Fonction d'estimation finale
  
  # Elle est la meme que dans ton code
  
  # ----------------------------------------------------------
  
  
  
  estimer <- function(sel, label) {
    
    tryCatch({
      
      Xs <- X[, sel, drop = FALSE]
      
      
      
      bM <- coef(
        
        cv.glmnet(cbind(A, M, Xs), Y, alpha = 1),
        
        "lambda.min"
        
      )["M", ]
      
      
      
      r_ATE <- dr_est(
        
        oal(
          
          Xs,
          
          A,
          
          nm((sc$S1[sel] * sc$S2[sel] * sc$S3[sel])^(1 / 3))
          
        ),
        
        Xs,
        
        A,
        
        Y
        
      )
      
      
      
      r_ACDE <- dr_est(
        
        oal(
          
          Xs,
          
          A,
          
          nm(sc$S3[sel])
          
        ),
        
        Xs,
        
        A,
        
        Y - bM * M
        
      )
      
      
      
      AIE <- r_ATE["est"] - r_ACDE["est"]
      
      se  <- sqrt(r_ATE["se"]^2 + r_ACDE["se"]^2)
      
      
      
      est <- c(r_ATE["est"], r_ACDE["est"], AIE)
      
      se2 <- c(r_ATE["se"],  r_ACDE["se"],  se)
      
      
      
      data.frame(
        
        effet    = c("ATE", "ACDE", "AIE"),
        
        estimate = est,
        
        se       = se2,
        
        lower    = est - z * se2,
        
        upper    = est + z * se2
        
      )
      
    }, error = function(e) {
      
      cat("ERREUR dans", label, ":", e$message, "\n")
      
      NULL
      
    })
    
  }
  
  
  
  list(
    
    deux_etapes = estimer(sel_2step, "deux_etapes"),
    
    geo = estimer(sel_geo, "geo"),
    
    maxZ = estimer(sel_maxZ, "maxZ"),
    
    
    
    hyb20_10_deux = estimer(sel_hyb20_10_deux, "hyb20_10_deux"),
    
    hyb20_10_maxZ = estimer(sel_hyb20_10_maxZ, "hyb20_10_maxZ"),
    
    
    
    hyb15_10_5_deux_S2 = estimer(sel_hyb15_10_5_deux_S2, "hyb15_10_5_deux_S2"),
    
    hyb15_10_5_maxZ_S2 = estimer(sel_hyb15_10_5_maxZ_S2, "hyb15_10_5_maxZ_S2")
    
  )
  
}