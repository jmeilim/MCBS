library(parallel)
library(Ball)
library(ggplot2)

# Charger les fonctions
source("MCBS/MCBS.R")
source("MCBS/DGP.R")

run_one <- function(b, n, p, lambda) {
  
  set.seed(1000 + b)
  
  dat <- DGP_rich(
    n = n,
    p = p,
    lambda = lambda
  )
  
  k <- floor(n / log(n))
  
  selected_geo <- screening(
    X = dat$X,
    Y = dat$Y,
    M = dat$M,
    A = dat$A,
    k = k
  )
  
  capture_all_confounders <- as.numeric(
    all(dat$true_confounders %in% selected_geo)
  )
  
  capture_pure_predictors <- mean(
    dat$pure_predictors %in% selected_geo
  )
  
  data.frame(
    capture_all_confounders = capture_all_confounders,
    capture_pure_predictors = capture_pure_predictors
  )
}

# -----------------------------
# Simulation parallèle
# -----------------------------
run_simulation <- function(
    B = 100,
    n,
    p,
    lambda,
    ncores = 16
) {
  
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, {
    library(Ball)
    source("MCBS/MCBS.R")
    NULL
  })
  
  clusterExport(
    cl,
    c("expit", "DGP_rich", "run_one", "n", "p", "lambda"),
    envir = environment()
  )
  
  all_res <- parLapply(cl, 1:B, function(b) {
    run_one(
      b = b,
      n = n,
      p = p,
      lambda = lambda
    )
  })
  
  all_res <- do.call(rbind, all_res)
  
  data.frame(
    prob_all_confounders = mean(all_res$capture_all_confounders),
    prob_pure_predictors = mean(all_res$capture_pure_predictors)
  )
}

# -----------------------------
# Paramètres
# -----------------------------
scenarios <- data.frame(
  n = c(300, 300, 600, 400, 500, 600),
  p = c(100, 1000, 200, 1000, 1500, 2000)
)

lambda_grid <- seq(0.5, 4, by = 0.5)

B <- 100
ncores <- 100

# -----------------------------
# Boucle principale
# -----------------------------
results <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  n_i <- scenarios$n[i]
  p_i <- scenarios$p[i]
  
  cat("\n=============================\n")
  cat("Scenario:", paste0("(", n_i, ", ", p_i, ")"), "\n")
  cat("=============================\n")
  
  for (lambda in lambda_grid) {
    
    cat("lambda =", lambda, "\n")
    
    res <- run_simulation(
      B = B,
      n = n_i,
      p = p_i,
      lambda = lambda,
      ncores = ncores
    )
    
    results <- rbind(
      results,
      data.frame(
        n = n_i,
        p = p_i,
        scenario = paste0("(", n_i, ", ", p_i, ")"),
        lambda = lambda,
        prob_all_confounders = res$prob_all_confounders,
        prob_pure_predictors = res$prob_pure_predictors
      )
    )
    
    print(tail(results, 1))
  }
}

write.csv(results, "results_rich_DGP_geo.csv", row.names = FALSE)

# -----------------------------
# Graphe principal
# -----------------------------
g <- ggplot(results, aes(
  x = lambda,
  y = prob_all_confounders,
  color = scenario
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = 0.8,
    color = "darkgreen",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Probabilité de récupérer tous les vrais confounders",
    subtitle = "Méthode geo — DGP enrichi avec signaux weak / medium / strong",
    x = expression(lambda),
    y = "P(tous les confounders sélectionnés)",
    color = "Scénario (n, p)"
  ) +
  theme_minimal(base_size = 13)

print(g)

ggsave(
  "graph_rich_DGP_geo.png",
  plot = g,
  width = 10,
  height = 6,
  dpi = 300
)