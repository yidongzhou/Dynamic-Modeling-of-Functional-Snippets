source("./code/gdm.R")
sapply(paste0('./mcfda/', list.files('mcfda')), source)
library(doSNOW)
cl <- makeCluster(50)
registerDoSNOW(cl)
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

gbmiter <- function(theta, sigma, delta, x0, z1){# function to generate sample paths of GBM
  x1 <- x0 * exp(delta * (theta - sigma^2 / 2) + sigma * sqrt(delta) * z1)
  x1
}

# GBM: dXt = theta Xt dt + sigma Xt dBt, t \in [0, 1]
theta <- 0.3
sigma <- 0.5
tau <- 1# length of the time interval
delta <- 0.05# time spacing
tgrid <- seq(0, 1, by = delta)# time grid
K <- tau / delta + 1# number of time points
x0 <- 1# initial value
mx <- x0 * exp(theta * tgrid)
sx <- x0 * exp(theta * tgrid) * sqrt(exp(sigma^2 * tgrid) - 1)
M <- 1000# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
Q <- 500# number of simulations
nvec <- c(50, 200, 1000)# number of subjects
nuvec <- c(0, 0.01, 0.1)# noise level

set.seed(1)
gbm1 <- list()
gbm2 <- list()
gbm3 <- list()
gbm4 <- list()
for(i in 1:length(nuvec)){
  gbm1[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  gbm2[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  gbm3[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  gbm4[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    gbmq <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
      x <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n){
        x[i, 1] <- x0
        for(j in 2:K){
          x[i, j] <- gbmiter(theta, sigma, delta, x[i, j - 1], rnorm(1))
        }
      }
      # add noises to the snippets
      x <- x + matrix(rnorm(n * K, sd = nu), nrow = n)
      
      Lt <- list()
      Ly <- list()
      for(i in 1:n){
        idx0 <- sample(K - 1, 1)
        Lt[[i]] <- tgrid[c(idx0, idx0 + 1)]
        Ly[[i]] <- x[i, c(idx0, idx0 + 1)]
      }
      
      z0 <- c(x0, 0)# (x0, t0)
      bm <- matrix(rnorm(M * (K - 1)), nrow = M)
      xtilde <- matrix(0, nrow = M, ncol = K)
      for(m in 1:M){
        xtilde[m, 1] <-  z0[1]
        for(j in 2:K){
          xtilde[m, j] <- gbmiter(theta, sigma, delta, xtilde[m, j - 1], bm[m, j - 1])
        }
      }
      xhat <- gdm(Ly = Ly, Lt = Lt, z0 = z0, tp = tgrid, optns = list(M = M, bm = bm))$path
      fitsp <- sp(Ly, Lt, bw = 0.1, domain = c(0, 1), newt = tgrid, x0, M, bm)
      c(sqrt(mean((xhat[, K] - xtilde[, K])^2)),# RMSE
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, mean) - mx)^2),# mean
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, sd) - sx)^2),# sd
        sqrt(mean((fitsp$path[, K] - xtilde[, K])^2)))# RMSE using Lin and Wang (2022) JASA
    }
    gbm1[[nui]][, ni] <- gbmq[, 1]
    gbm2[[nui]][, ni] <- gbmq[, 2]
    gbm3[[nui]][, ni] <- gbmq[, 3]
    gbm4[[nui]][, ni] <- gbmq[, 4]
  }
}
stopCluster(cl)
save(gbm1, gbm2, gbm3, gbm4, file = "./data/gbm.RData")
# sapply(gbm1, colMeans)# RMSE
# sapply(gbm2, colMeans)# mean
# sapply(gbm3, colMeans)# sd
# sapply(gbm4, colMeans)# RMSE using Lin and Wang (2022) JASA
