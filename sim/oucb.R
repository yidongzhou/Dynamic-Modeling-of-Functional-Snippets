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
nvec <- c(50, 200, 1000)# number of subjects c(100, 200, 500, 1000, 2000, 5000)
nuvec <- c(0, 0.01, 0.1)# noise level
alpha <- 0.05 # confidence level

ou <- list()
for(i in 1:length(nuvec)){
  ou[[i]] <- matrix(0, nrow = Q, ncol = length(nvec))
}
for(nui in 1:length(nuvec)){
  nu <- nuvec[nui]
  for(ni in 1:length(nvec)) {
    n <- nvec[ni]
    ouq <- foreach(1:Q, .combine = c, .options.snow = list(progress = progress)) %dopar% {
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
      xhat <- gdm(Ly = Ly, Lt = Lt, z0 = z0, tp = tgrid, optns = list(M = 2 * M))$path
      ci <- quantile(xhat[1:M, K], probs =c(alpha / 2, (1 - alpha / 2)))
      sum(xhat[(M + 1):(2 * M), K] >= ci[1] & xhat[(M + 1):(2 * M), K] <= ci[2]) / M# coverage
    }
    ou[[nui]][, ni] <- ouq
  }
}
stopCluster(cl)
save(ou, file = "data/oucb.RData")