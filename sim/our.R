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

ouiter <- function(theta, sigma, delta, x0, z1){# function to generate sample paths of OU
  x1 <- x0 * exp(-theta * delta) + 
    sqrt(sigma^2 * (1 - exp(-2 * theta * delta)) / (2 * theta)) * z1
  x1
}

# OU process: dXt = - theta Xt dt + sigma dBt, t \in [0, 1]
theta <- 1
sigma <- 1
tau <- 1# length of the time interval
delta <- 0.05# time spacing
tgrid <- seq(0, 1, by = delta)# time grid
K <- tau / delta + 1# number of time points
x0 <- 0# initial value
mx <- x0 * exp(-theta * tgrid)
sx <- sqrt(sigma^2*(1 - exp(-2 * theta * tgrid)) / (2 * theta))
M <- 1000# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
Q <- 500# number of simulations
nvec <- c(50, 200, 1000)# number of subjects

#########################
### Five measurements ###
#########################
set.seed(1)
ou <- list()
for(i in 1:3){
  ou[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(ni in 1:length(nvec)) {
  n <- nvec[ni]
  ouq <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
    x <- matrix(0, nrow = n, ncol = K)
    for(i in 1:n){
      x[i, 1] <- x0
      for(j in 2:K){
        x[i, j] <- ouiter(theta, sigma, delta, x[i, j - 1], rnorm(1))
      }
    }
    
    Lt <- list()
    Ly <- list()
    Lt0 <- list()
    Ly0 <- list()
    for(i in 1:n){
      idx0 <- sample(K - 4, 1)
      Lt[[i]] <- tgrid[idx0:(idx0 + 4)]
      Ly[[i]] <- x[i, idx0:(idx0 + 4)]
      idx <- sample(4, 1)# randomly select two consecutive observations
      Lt0[[i]] <- Lt[[i]][c(idx, idx + 1)]
      Ly0[[i]] <- Ly[[i]][c(idx, idx + 1)]
    }
    
    z0 <- c(x0, 0)# (x0, t0)
    bm <- matrix(rnorm(M * (K - 1)), nrow = M)
    xtilde <- matrix(0, nrow = M, ncol = K)
    for(m in 1:M){
      xtilde[m, 1] <-  z0[1]
      for(j in 2:K){
        xtilde[m, j] <- ouiter(theta, sigma, delta, xtilde[m, j - 1], bm[m, j - 1])
      }
    }
    xhat <- gdm(Ly = Ly, Lt = Lt, z0 = z0, tp = tgrid, optns = list(M = M, bm = bm))$path
    xhat0 <- gdm(Ly = Ly0, Lt = Lt0, z0 = z0, tp = tgrid, optns = list(M = M, bm = bm))$path
    fitsp <- sp(Ly, Lt, bw = 0.1, domain = c(0, 1), newt = tgrid, x0, M, bm)
    c(sqrt(mean((xhat[, K] - xtilde[, K])^2)),# RMSE
      sqrt(mean((xhat0[, K] - xtilde[, K])^2)),# RMSE using only two measurements
      sqrt(mean((fitsp$path[, K] - xtilde[, K])^2)))# RMSE using Lin and Wang (2022) JASA
  }
  ou[[1]][, ni] <- ouq[, 1]
  ou[[2]][, ni] <- ouq[, 2]
  ou[[3]][, ni] <- ouq[, 3]
}
stopCluster(cl)
save(ou, file = "./data/our.RData")
# sapply(ou, colMeans)
