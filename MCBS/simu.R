library(parallel)
library(Ball)
library(ggplot2)

# Charger les fonctions
source("MCBS/MCBS.R")
source("MCBS/DGP.R")

# -----------------------------
# Une simulation : geo sélectionne-t-il X2 ?
# -----------------------------
run_one_x2_geo <- function(b, n, p, k, beta_AY, beta_MY, beta_AM) {
  
  set.seed(1000 + b)
  
  dat <- DGP(
    n = n,
    p = p,
    beta_AY = beta_AY,
    beta_MY = beta_MY,
    beta_AM = beta_AM
  )
  
  selected_geo <- screening(
    X = dat$X,
    Y = dat$Y,
    M = dat$M,
    A = dat$A,
    k = k
  )
  
  # X2 = variable 2
  as.numeric(2 %in% selected_geo)
}

# -----------------------------
# Pour un scénario (n,p) et une valeur beta_AM
# retourne la proba de sélectionner X2
# -----------------------------
run_prob_x2_geo <- function(
    B = 100,
    n,
    p,
    k,
    beta_AY,
    beta_MY,
    beta_AM,
    ncores = 100
) {
  
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, {
    library(Ball)
    source("MCBS/MCBS.R")
    source("MCBS/DGP.R")
    NULL
  })
  
  clusterExport(
    cl,
    c("run_one_x2_geo", "n", "p", "k", "beta_AY", "beta_MY", "beta_AM"),
    envir = environment()
  )
  
  sel_x2 <- unlist(
    parLapply(cl, 1:B, function(b) {
      run_one_x2_geo(
        b = b,
        n = n,
        p = p,
        k = k,
        beta_AY = beta_AY,
        beta_MY = beta_MY,
        beta_AM = beta_AM
      )
    })
  )
  
  mean(sel_x2)
}

# -----------------------------
# Grille des scénarios
# -----------------------------
scenarios <- data.frame(
  n = c(300, 300, 600, 400, 500, 600),
  p = c(100, 1000, 200, 1000, 1500, 2000)
)

# -----------------------------
# Grille des beta_AM
# -----------------------------
beta_grid <- seq(0.1, 1.5, by = 0.1)

# -----------------------------
# Paramètres fixes
# -----------------------------
B <- 100
k <- 30
beta_AY <- 0.2
beta_MY <- 0.6
ncores <- 16   # adapte selon ton cluster

# -----------------------------
# Boucle principale
# -----------------------------
results <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  n_i <- scenarios$n[i]
  p_i <- scenarios$p[i]
  
  cat("\n=============================\n")
  cat("Scenario:", "(n,p) =", paste0("(", n_i, ",", p_i, ")"), "\n")
  cat("=============================\n")
  
  for (beta in beta_grid) {
    
    cat("beta_AM =", beta, "\n")
    
    prob_x2 <- run_prob_x2_geo(
      B = B,
      n = n_i,
      p = p_i,
      k = k,
      beta_AY = beta_AY,
      beta_MY = beta_MY,
      beta_AM = beta,
      ncores = ncores
    )
    
    results <- rbind(
      results,
      data.frame(
        n = n_i,
        p = p_i,
        scenario = paste0("(", n_i, ", ", p_i, ")"),
        beta_AM = beta,
        prob_x2 = prob_x2
      )
    )
    
    cat("   -> P(selection X2) =", prob_x2, "\n")
  }
}

print(results)

# Sauvegarde éventuelle
write.csv(results, "results_prob_X2_geo.csv", row.names = FALSE)

# -----------------------------
# Seuil beta_AM minimal pour atteindre 0.8
# -----------------------------
thresholds <- do.call(rbind, lapply(split(results, results$scenario), function(df) {
  
  df <- df[order(df$beta_AM), ]
  idx <- which(df$prob_x2 >= 0.8)
  
  if (length(idx) == 0) {
    beta_threshold <- NA
  } else {
    beta_threshold <- df$beta_AM[min(idx)]
  }
  
  data.frame(
    scenario = df$scenario[1],
    beta_threshold_0.8 = beta_threshold
  )
}))

print(thresholds)

write.csv(thresholds, "thresholds_X2_geo_0.8.csv", row.names = FALSE)

# -----------------------------
# Graphe
# -----------------------------
g <- ggplot(results, aes(x = beta_AM, y = prob_x2, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = beta_grid) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Probabilité de sélectionner X2 selon beta_AM",
    subtitle = "Méthode geo — une courbe par scénario (n, p)",
    x = expression(beta[AM]),
    y = "P(X2 sélectionné)",
    color = "Scénario (n, p)"
  ) +
  theme_minimal(base_size = 13)

print(g)

ggsave("graph_prob_X2_geo.png", plot = g, width = 10, height = 6, dpi = 300)