source("MCBS/MCBS.R")
source("MCBS/DGP.R")  

test_rich <- function(n=300, p=1000, lambda=1.5, B=50, k=30) {
  
  
  CAY_2step <- CAM_2step <- CMY_2step <- PY_2step <- PM_2step <- numeric(B)
  CAY_geo   <- CAM_geo   <- CMY_geo  <-  PY_geo <-  PM_geo   <- numeric(B)
  
  for(b in 1:B) {
    set.seed(1000 + b)
    dat <- DGP_rich(n, p, lambda)
    
    s2 <- screening_deux_etapes(dat$X, dat$Y, dat$M, dat$A, k)
    sg <- screening_geo(dat$X, dat$Y, dat$M, dat$A, k)
    
    # C_AY = indices 1,2,3
    # C_AM = indices 4,5,6
    # C_MY = indices 7,8,9
    CAY_2step[b] <- mean(1:3 %in% s2)
    CAM_2step[b] <- mean(4:6 %in% s2)
    CMY_2step[b] <- mean(7:9 %in% s2)
    PY_2step[b] <- mean(10:12 %in% s2)
    PM_2step[b] <- mean(13:15 %in% s2)
    CAY_geo[b]   <- mean(1:3 %in% sg)
    CAM_geo[b]   <- mean(4:6 %in% sg)
    CMY_geo[b]   <- mean(7:9 %in% sg)
    PY_geo[b] <- mean(10:12 %in% sg)
    PM_geo[b] <- mean(13:15 %in% sg)
  }
  
  cat("\n=== DGP_rich n=",n,"p=",p,"lambda=",lambda,"===\n")
  cat("         2step    geo\n")
  cat("C_AY  :", round(mean(CAY_2step),3),
      round(mean(CAY_geo),3), "\n")
  cat("C_AM  :", round(mean(CAM_2step),3),
      round(mean(CAM_geo),3), "\n")
  cat("C_MY  :", round(mean(CMY_2step),3),
      round(mean(CMY_geo),3), "\n") 
  cat("PY  :", round(mean(PY_2step),3),
      round(mean(PY_geo),3), "\n")
  cat("PM  :", round(mean(PY_2step),3),
      round(mean(PM_geo),3), "\n")
}


t0 <- Sys.time()
test_rich(n=300, p=1000, lambda=4, B=50)
test_rich(n=500, p=1000, lambda=4, B=50)
cat("Temps:", round(Sys.time()-t0,1), "sec\n")