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

hliter <- function(sigma, t0, t1, x0, z1){# function to generate sample paths of HL
  x1 <- x0 + sin(t1) - sin(t0) + sigma * sqrt(t1 - t0) * z1
  x1
}

# HL process: dXt = cos(t) dt + sigma dBt, t \in [0, 1]
sigma <- 1
tau <- 1# length of the time interval
delta <- 0.05# time spacing
tgrid <- seq(0, 1, by = delta)# time grid
K <- tau / delta + 1# number of time points
x0 <- 0# initial value
mx <- x0 + sin(tgrid)
sx <- sigma * sqrt(tgrid)
M <- 1000# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
Q <- 500# number of simulations
nvec <- c(50, 200, 1000)# number of subjects
nuvec <- c(0, 0.01, 0.1)# noise level

set.seed(1)
hl1 <- list()
hl2 <- list()
hl3 <- list()
hl4 <- list()
for(i in 1:length(nuvec)){
  hl1[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  hl2[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  hl3[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
  hl4[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    hlq <- foreach(1:Q, .combine = rbind, .options.snow = list(progress = progress)) %dopar% {
      x <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n){
        x[i, 1] <- x0
        for(j in 2:K){
          x[i, j] <- hliter(sigma, tgrid[j - 1], tgrid[j], x[i, j - 1], rnorm(1))
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
          xtilde[m, j] <- hliter(sigma, tgrid[j - 1], tgrid[j], xtilde[m, j - 1], bm[m, j - 1])
        }
      }
      xhat <- gdm(Ly = Ly, Lt = Lt, z0 = z0, tp = tgrid, optns = list(M = M, bm = bm))$path
      fitsp <- sp(Ly, Lt, bw = 0.1, domain = c(0, 1), newt = tgrid, x0, M, bm)
      c(sqrt(mean((xhat[, K] - xtilde[, K])^2)),# RMSE
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, mean) - mx)^2),# mean
        fdapace::trapzRcpp(tgrid, (apply(xhat, 2, sd) - sx)^2),# sd
        sqrt(mean((fitsp$path[, K] - xtilde[, K])^2)))# RMSE using Lin and Wang (2022) JASA
    }
    hl1[[nui]][, ni] <- hlq[, 1]
    hl2[[nui]][, ni] <- hlq[, 2]
    hl3[[nui]][, ni] <- hlq[, 3]
    hl4[[nui]][, ni] <- hlq[, 4]
  }
}
stopCluster(cl)
save(hl1, hl2, hl3, hl4, file = "./data/hl.RData")
# sapply(hl1, colMeans)# RMSE
# sapply(hl2, colMeans)# mean
# sapply(hl3, colMeans)# sd
# sapply(hl4, colMeans)# RMSE using Lin and Wang (2022) JASA
