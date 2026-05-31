library(causalweight)
library(MASS)
set.seed(123)

n <- 1000
p <- 50

Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X) <- paste0("X", 1:p)

beta <- 0.3 / (1:p)^2

D <- as.numeric(X %*% beta + rnorm(n) > 0)
M <- as.numeric(0.5 * D + X %*% beta + rnorm(n) > 0)
Y <- as.numeric(0.5 * D + 0.5 * M + 0.5 * D * M + X %*% beta + rnorm(n))

res <- medDML(y = Y, d = D, m = M, x = X, trim = 0.05)
print(res)