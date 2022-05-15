library("devtools")
library("dplyr")
library(nloptr)
help(nloptr)
N <- 300
K <- 20
lb <- rep(1, N)
ub <- rep(3, N)
objective <- function(x, K0) {
    return(-1 * sum(tail(x, K0)) / sum(x))
}
x0 <- rep(2, N)


res <- slsqp(x0, objective, lower=lb, upper=ub, K0=K)
print(res$par)
print(-res$value)
print(-objective(c(rep(1, N - K), rep(3, K)), K))

#%%-----------------------------------------------------------------------


set.seed(123)
N <- 300
K <- 5
lb <- rep(1, N)
ub <- rep(3, N)
objective <- function(x, K0) {
    return(-1 * sum(tail(x, K0)) / sum(x))
}

coef <- rnorm(N, 0, 0.2)
con <- function(x, coef0=coef) {
    return(sum(x * coef0) - 3)
}

x0 <- runif(N, 1, 3)
x0 <- rep(1, N)
x0 <- c(rep(1, N - K), rep(3, K))
res <- slsqp(x0, objective, heq = con, lower = lb, upper = ub, K0 = K, 
        control = list(maxeval = 1000, xtol_rel = 1e-5, ftol_abs = 1e-3))

print(res$par)
print(-res$value)
print(-objective(c(rep(1, N - K), rep(3, K)), K))
print(con(res$par, coef))
print(res$message)
print(res$iter)




