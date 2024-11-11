library(doSNOW)
cl <- makeCluster(8)
registerDoSNOW(cl)
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

M <- 1000# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
Q <- 500# number of simulations
nvec <- c(50, 200, 1000)# number of subjects
nuvec <- c(0, 0.01, 0.1)# noise level
tau <- 1# length of the time interval
delta <- 0.05# time spacing
tgrid <- seq(0, 1, by = delta)# time grid
K <- tau / delta + 1# number of time points

# HL process: dXt = cos(t) dt + sigma dBt, t \in [0, 1]
hliter <- function(sigma, t0, t1, x0, z1){# function to generate sample paths of HL
  x1 <- x0 + sin(t1) - sin(t0) + sigma * sqrt(t1 - t0) * z1
  x1
}
sigma <- 1
x0 <- 0# initial value
mx <- x0 + sin(tgrid)
sx <- sigma * sqrt(tgrid)

set.seed(1)
hl <- list()
for(i in 1:length(nuvec)){
  hl[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    hl[[nui]][, ni] <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
      x <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n){
        x[i, 1] <- x0
        for(j in 2:K){
          x[i, j] <- hliter(sigma, tgrid[j - 1], tgrid[j], x[i, j - 1], rnorm(1))
        }
      }
      # add noises to the snippets
      x <- x + matrix(rnorm(n * K, sd = nu), nrow = n)
      
      df <- matrix(0, nrow = n * 2, ncol = 2)
      for(i in 1:n){
        idx0 <- sample(K - 1, 1)
        df[(2 * (i - 1) + 1) : (2 * i), 1] <- tgrid[c(idx0, idx0 + 1)]
        df[(2 * (i - 1) + 1) : (2 * i), 2] <- x[i, c(idx0, idx0 + 1)]
      }
      df <- df[order(df[, 1]), ]
      y <- NULL
      for(tk in unique(df[, 1])) {
        y <- c(y, mean(df[df[, 1] == tk, 2]))
      }
      
      fit <- forecast::auto.arima(y)
      # fit <- arima(y, order = c(1, 0, 1))
      # psi1 <- coef(fit)[1]
      # theta1 <- coef(fit)[2]
      # alpha <- coef(fit)[3]
      # sigma2 <- fit$sigma2
      
      bm <- matrix(rnorm(M * K), nrow = M)
      xtilde <- matrix(0, nrow = M, ncol = K)
      xhat <- matrix(0, nrow = M, ncol = K)
      for(m in 1:M){
        xhat[m, ] <- simulate(fit, future = FALSE, innov = bm[m, ])
        xtilde[m, 1] <-  x0
        for(j in 2:K){
          xtilde[m, j] <- hliter(sigma, tgrid[j - 1], tgrid[j], xtilde[m, j - 1], bm[m, j - 1])
          # xhat[m, j] <- alpha + psi1 * xtilde[m, j - 1] + theta1 * sqrt(sigma2) * bm[m, j - 1] + sqrt(sigma2) * bm[m, j]
        }
      }
      sqrt(mean((xhat[, K] - xtilde[, K])^2))# RMSE
    }
  }
}

# OU process: dXt = - theta Xt dt + sigma dBt, t \in [0, 1]
ouiter <- function(theta, sigma, delta, x0, z1){# function to generate sample paths of OU
  x1 <- x0 * exp(-theta * delta) + 
    sqrt(sigma^2 * (1 - exp(-2 * theta * delta)) / (2 * theta)) * z1
  x1
}
theta <- 1
sigma <- 1
x0 <- 0# initial value
mx <- x0 * exp(-theta * tgrid)
sx <- sqrt(sigma^2 * (1 - exp(-2 * theta * tgrid)) / (2 * theta))

set.seed(1)
ou <- list()
for(i in 1:length(nuvec)){
  ou[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    ou[[nui]][, ni] <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
      x <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n){
        x[i, 1] <- x0
        for(j in 2:K){
          x[i, j] <- ouiter(theta, sigma, delta, x[i, j - 1], rnorm(1))
        }
      }
      # add noises to the snippets
      x <- x + matrix(rnorm(n * K, sd = nu), nrow = n)
      
      df <- matrix(0, nrow = n * 2, ncol = 2)
      for(i in 1:n){
        idx0 <- sample(K - 1, 1)
        df[(2 * (i - 1) + 1) : (2 * i), 1] <- tgrid[c(idx0, idx0 + 1)]
        df[(2 * (i - 1) + 1) : (2 * i), 2] <- x[i, c(idx0, idx0 + 1)]
      }
      df <- df[order(df[, 1]), ]
      y <- NULL
      for(tk in unique(df[, 1])) {
        y <- c(y, mean(df[df[, 1] == tk, 2]))
      }
      
      fit <- forecast::auto.arima(y)
      # fit <- arima(y, order = c(1, 1, 1))
      # psi1 <- coef(fit)[1]
      # theta1 <- coef(fit)[2]
      # alpha <- coef(fit)[3]
      # sigma2 <- fit$sigma2
      
      bm <- matrix(rnorm(M * K), nrow = M)
      xtilde <- matrix(0, nrow = M, ncol = K)
      xhat <- matrix(0, nrow = M, ncol = K)
      for(m in 1:M){
        xhat[m, ] <- simulate(fit, future = FALSE, innov = bm[m, ])
        xtilde[m, 1] <-  x0
        for(j in 2:K){
          xtilde[m, j] <- ouiter(theta, sigma, delta, xtilde[m, j - 1], bm[m, j - 1])
          # xhat[m, j] <- alpha + psi1 * xtilde[m, j - 1] + theta1 * sqrt(sigma2) * bm[m, j - 1] + sqrt(sigma2) * bm[m, j]
        }
      }
      sqrt(mean((xhat[, K] - xtilde[, K])^2))# RMSE
    }
  }
}

# GBM: dXt = theta Xt dt + sigma Xt dBt, t \in [0, 1]
gbmiter <- function(theta, sigma, delta, x0, z1){# function to generate sample paths of GBM
  x1 <- x0 * exp(delta * (theta - sigma^2 / 2) + sigma * sqrt(delta) * z1)
  x1
}
theta <- 0.3
sigma <- 0.5
x0 <- 1# initial value
mx <- x0 * exp(theta * tgrid)
sx <- x0 * exp(theta * tgrid) * sqrt(exp(sigma^2 * tgrid) - 1)

set.seed(1)
gbm <- list()
for(i in 1:length(nuvec)){
  gbm[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    gbm[[nui]][, ni] <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
      x <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n){
        x[i, 1] <- x0
        for(j in 2:K){
          x[i, j] <- gbmiter(theta, sigma, delta, x[i, j - 1], rnorm(1))
        }
      }
      # add noises to the snippets
      x <- x + matrix(rnorm(n * K, sd = nu), nrow = n)
      
      df <- matrix(0, nrow = n * 2, ncol = 2)
      for(i in 1:n){
        idx0 <- sample(K - 1, 1)
        df[(2 * (i - 1) + 1) : (2 * i), 1] <- tgrid[c(idx0, idx0 + 1)]
        df[(2 * (i - 1) + 1) : (2 * i), 2] <- x[i, c(idx0, idx0 + 1)]
      }
      df <- df[order(df[, 1]), ]
      y <- NULL
      for(tk in unique(df[, 1])) {
        y <- c(y, mean(df[df[, 1] == tk, 2]))
      }
      
      fit <- forecast::auto.arima(y)
      # fit <- arima(y, order = c(1, 0, 1))
      # psi1 <- coef(fit)[1]
      # theta1 <- coef(fit)[2]
      # alpha <- coef(fit)[3]
      # sigma2 <- fit$sigma2
      
      bm <- matrix(rnorm(M * K), nrow = M)
      xtilde <- matrix(0, nrow = M, ncol = K)
      xhat <- matrix(0, nrow = M, ncol = K)
      for(m in 1:M){
        xhat[m, ] <- simulate(fit, future = FALSE, innov = bm[m, ])
        xtilde[m, 1] <-  x0
        for(j in 2:K){
          xtilde[m, j] <- gbmiter(theta, sigma, delta, xtilde[m, j - 1], bm[m, j - 1])
          # xhat[m, j] <- alpha + psi1 * xtilde[m, j - 1] + theta1 * sqrt(sigma2) * bm[m, j - 1] + sqrt(sigma2) * bm[m, j]
        }
      }
      sqrt(mean((xhat[, K] - xtilde[, K])^2))# RMSE
    }
  }
}
stopCluster(cl)
save(hl, ou, gbm, file = "data/arma.RData")
# sapply(ou, colMeans)# OU
# sapply(hl, colMeans)# HL
# sapply(gbm, colMeans)# GBM