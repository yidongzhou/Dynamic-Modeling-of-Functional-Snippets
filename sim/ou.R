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
sx <- sqrt(sigma^2 * (1 - exp(-2 * theta * tgrid)) / (2 * theta))
M <- 1000# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
Q <- 500# number of simulations
nvec <- c(50, 200, 1000)# number of subjects
nuvec <- c(0, 0.01, 0.1)# noise level

set.seed(1)
ou1 <- list()
ou2 <- list()
ou3 <- list()
ou4 <- list()
for(i in 1:length(nuvec)){
  ou1[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  ou2[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  ou3[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  ou4[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
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
          xtilde[m, j] <- ouiter(theta, sigma, delta, xtilde[m, j - 1], bm[m, j - 1])
        }
      }
      xhat <- gdm(Ly = Ly, Lt = Lt, z0 = z0, tp = tgrid, optns = list(M = M, bm = bm))$path
      fitsp <- sp(Ly, Lt, bw = 0.1, domain = c(0, 1), newt = tgrid, x0, M, bm)
      c(sqrt(mean((xhat[, K] - xtilde[, K])^2)),# RMSE
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, mean) - mx)^2),# mean
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, sd) - sx)^2),# sd
        sqrt(mean((fitsp$path[, K] - xtilde[, K])^2)))# RMSE using Lin and Wang (2022) JASA
    }
    ou1[[nui]][, ni] <- ouq[, 1]
    ou2[[nui]][, ni] <- ouq[, 2]
    ou3[[nui]][, ni] <- ouq[, 3]
    ou4[[nui]][, ni] <- ouq[, 4]
  }
}
stopCluster(cl)
save(ou1, ou2, ou3, ou4, file = "./data/ou.RData")
# sapply(ou1, colMeans)# RMSE
# sapply(ou2, colMeans)# mean
# sapply(ou3, colMeans)# sd
# sapply(ou4, colMeans)# RMSE using Lin and Wang (2022) JASA
