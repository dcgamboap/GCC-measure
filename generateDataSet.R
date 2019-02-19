
library(MASS)

varianceMatrix1 <- matrix(0.9, nrow = 100, ncol = 100)
diag(varianceMatrix1) <- 1


varianceMatrix2 <- matrix(0.9, nrow = 100, ncol = 100)
diag(varianceMatrix2) <- 1

x1 <- mvrnorm(n = 1000, mu = rep(0, 100), Sigma = varianceMatrix1)
x2 <- mvrnorm(n = 1000, mu = rep(0, 100), Sigma = varianceMatrix2)
x3 <- mvrnorm(n = 1000, mu = rep(0, 30), Sigma = diag(1, nrow = 30))


x = cbind(x1, x2, x3)
