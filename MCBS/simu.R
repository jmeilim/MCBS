source("MCBS/MCBS.R")

set.seed(42)

run_one <- function(n, p, k, beta_AY, beta_MY, beta_AM) {
  
  data <- DGP(n, p, beta_AY = beta_AY, beta_MY = beta_MY, beta_AM = beta_AM)
  
  invisible(capture.output({
    sel_geo <- screening(data$X, data$Y, data$M, data$A, k)
    sel_s1  <- screening_S1_clean(data$X, data$Y, data$A, k)
  }))
  
  res <- data.frame(
    method = c("geo", "geo", "geo", "S1", "S1", "S1"),
    confounder = c("X1_C_AY", "X2_C_AM", "X3_C_MY",
                   "X1_C_AY", "X2_C_AM", "X3_C_MY"),
    selected = c(
      1 %in% sel_geo,
      2 %in% sel_geo,
      3 %in% sel_geo,
      1 %in% sel_s1,
      2 %in% sel_s1,
      3 %in% sel_s1
    )
  )
  
  return(res)
}
run_simulation <- function(B = 100, n, p, k, beta_AY, beta_MY, beta_AM) {
  
  all_res <- list()
  
  for(b in 1:B) {
    
    all_res[[b]] <- run_one(
      n = n,
      p = p,
      k = k,
      beta_AY = beta_AY,
      beta_MY = beta_MY,
      beta_AM = beta_AM
    )
  }
  
  all_res <- do.call(rbind, all_res)
  
  proba <- aggregate(
    selected ~ method + confounder,
    data = all_res,
    FUN = mean
  )
  
  return(proba)
}
for(beta in c(0.05, 0.1, 0.2, 0.6, 0.9)) {
  cat("\n beta_AM =", beta, "\n")
  
  res <- run_simulation(
    B = 100,
    n = 200,
    p = 20,
    k = 5,
    beta_AY = 0.2,
    beta_MY = 0.6,
    beta_AM = beta
  )
  
  print(res)
}