## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(Matrix)
library(WoodburyMatrix)

set.seed(20200721)

A1 <- rsparsematrix(500, 500, 0.05)
B1 <- rsparsematrix(500, 500, 0.05)
W1 <- WoodburyMatrix(A1, B1)

## -----------------------------------------------------------------------------
A2 <- crossprod(rsparsematrix(500, 500, 0.05))
B2 <- crossprod(rsparsematrix(500, 500, 0.05))
W2 <- WoodburyMatrix(A2, B2, symmetric = TRUE)

## -----------------------------------------------------------------------------
n <- 2000
p <- 20
sigma_squared <- 1
X <- matrix(rnorm(n * p), nrow = n)
Q_beta <- Diagonal(p)
mu_beta <- rep(0, p)

Sigma_Y <- WoodburyMatrix(
  Diagonal(n, 1 / sigma_squared),
  Q_beta,
  X = X
)

## -----------------------------------------------------------------------------
Sigma_Y_direct <- instantiate(Sigma_Y)
object.size(Sigma_Y)
object.size(instantiate(Sigma_Y))

## -----------------------------------------------------------------------------
Y <- rwnorm(1, X %*% mu_beta, Sigma_Y)
system.time(print(dwnorm(Y, X %*% mu_beta, Sigma_Y, log = TRUE)))

## -----------------------------------------------------------------------------
Y_m <- Y - X %*% mu_beta
system.time(print(as.numeric(
  - 0.5 * n * log(2 * pi)
  - 0.5 * determinant(Sigma_Y_direct)$modulus
  - 0.5 * crossprod(Y_m, solve(Sigma_Y_direct, Y_m))
)))

## -----------------------------------------------------------------------------
n <- 10000
rho <- 0.95
Q_1 <- bandSparse(
  n,
  diagonals = list(
    c(1, rep(1 + rho ^ 2, n - 2), 1),
    rep(-rho, n - 1)
  ),
  k = c(0, 1),
  symmetric = TRUE
)
Q_2 <- Diagonal(n)

Sigma_Y <- WoodburyMatrix(Q_1, Q_2)

## -----------------------------------------------------------------------------
print(isSymmetric(Sigma_Y))

## -----------------------------------------------------------------------------
Y <- rwnorm(1, covariance = Sigma_Y)
plot(Y, type = 'l')

## -----------------------------------------------------------------------------
system.time(print(dwnorm(Y, covariance = Sigma_Y, log = TRUE)))

## -----------------------------------------------------------------------------
n <- 500
p1 <- 50
p2 <- 10
A <- rsparsematrix(n, n, 0.5)
B1 <- rsparsematrix(p1, p1, 0.9)
B2 <- rsparsematrix(p2, p2, 0.9)
U1 <- rsparsematrix(n, p1, 0.9)
U2 <- rsparsematrix(n, p2, 0.9)
W <- WoodburyMatrix(A, list(B1, B2), U = list(U1, U2), V = list(t(U1), t(U2)))
b <- rnorm(n)
str(solve(W, b))

