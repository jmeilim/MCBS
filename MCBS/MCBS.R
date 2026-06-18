# ============================================================
# MCBS.R — CBS Mediation
# ============================================================
library(Ball)
library(glmnet)

# ── Scores ──────────────────────────────────────────────────

S1 <- function(x, Y, A) {
  idx0 <- A==0; idx1 <- A==1
  mean(idx0)*bcov(x[idx0], Y[idx0])^2 +
    mean(idx1)*bcov(x[idx1], Y[idx1])^2
}

S2 <- function(x, M, A) S1(x, M, A)

S3 <- function(x, Y, A, M) {
  grps <- list(A==0&M==0, A==0&M==1, A==1&M==0, A==1&M==1)
  out  <- 0
  for(g in grps) {
    if(sum(g) < 2) next
    out <- out + mean(g)*bcov(x[g], Y[g])^2
  }
  out
}

# ── Screening ───────────────────────────────────────────────

compute_scores <- function(X, Y, M, A) {
  p <- ncol(X)
  s1 <- s2 <- s3 <- numeric(p)
  for(j in 1:p) {
    s1[j] <- S1(X[,j], Y, A)
    s2[j] <- S2(X[,j], M, A)
    s3[j] <- S3(X[,j], Y, A, M)
  }
  list(S1=s1, S2=s2, S3=s3)
}

screening_geo <- function(X, Y, M, A, k) {
  sc <- compute_scores(X, Y, M, A)
  order((rank(-sc$S1)*rank(-sc$S2)*rank(-sc$S3))^(1/3))[1:k]
}

screening_deux_etapes <- function(X, Y, M, A, k) {
  sc  <- compute_scores(X, Y, M, A)
  k2  <- floor(k/3); k1 <- k - k2
  top <- order(-sc$S2)[1:k2]
  rst <- setdiff(1:ncol(X), top)
  geo <- order((rank(-sc$S1[rst])*rank(-sc$S3[rst]))^0.5)[1:k1]
  union(top, rst[geo])
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

CBS_mediation <- function(X, A, M, Y, k=30, alpha=0.05) {
  p  <- ncol(X); k <- min(k,p)
  sc <- compute_scores(X, Y, M, A)
  
  # Screening deux étapes
  k2  <- floor(k/3); k1 <- k-k2
  top <- order(-sc$S2)[1:k2]
  rst <- setdiff(1:p, top)
  sel <- union(top,
               rst[order((rank(-sc$S1[rst])*rank(-sc$S3[rst]))^.5)[1:k1]])
  
  Xs <- X[,sel]
  nm <- function(s) s/(max(s)+1e-10)
  
  # OAL + DR
  bM    <- coef(cv.glmnet(cbind(A,M,Xs),Y,alpha=1),"lambda.min")["M",]
  r_ATE  <- dr_est(oal(Xs,A,nm((sc$S1[sel]*sc$S2[sel]*sc$S3[sel])^(1/3))),Xs,A,Y)
  r_ACDE <- dr_est(oal(Xs,A,nm(sc$S3[sel])), Xs,A,Y-bM*M)
  
  AIE <- r_ATE["est"] - r_ACDE["est"]
  se  <- sqrt(r_ATE["se"]^2 + r_ACDE["se"]^2)
  z   <- qnorm(1-alpha/2)
  
  data.frame(
    effet    = c("ATE","ACDE","AIE"),
    estimate = c(r_ATE["est"], r_ACDE["est"], AIE),
    se       = c(r_ATE["se"],  r_ACDE["se"],  se),
    lower    = c(r_ATE["est"], r_ACDE["est"], AIE) - z*c(r_ATE["se"],r_ACDE["se"],se),
    upper    = c(r_ATE["est"], r_ACDE["est"], AIE) + z*c(r_ATE["se"],r_ACDE["se"],se)
  )
}