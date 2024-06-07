library(ggplot2)
source('code/gdm.R')
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
M <- 100# number of Monte Carlo iterations to estimate the strong error E|Xhat- Xtilde|
n <- 200# number of subjects c(100, 200, 500, 1000, 2000, 5000)
nu <- 0.1# noise level

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

# load("data/ouplot.RData")
pdf(file = "latex/img/ou.pdf", width = 12, height = 8)
ggplot(data = data.frame(y = c(x, unlist(Ly), xtilde, xhat),
                         x = c(rep(tgrid, each = n), unlist(Lt), rep(rep(tgrid, each = M), 2)),
                         group = factor(c(rep(1:n, K), rep(1:n, each = 2), rep(rep(1:M, K), 2))),
                         class = factor(rep(c("simulated sample paths", "simulated snippets", 
                                              "true sample paths", "estimated sample paths"), 
                                            c(n * K, n * 2, M * K, M * K)),
                                        levels = c("simulated sample paths", "simulated snippets", 
                                                   "true sample paths", "estimated sample paths"))),
       aes(x= x, y = y, group = group, color = class)) +
  geom_line(data = ~subset(.x, class == "simulated snippets"), linewidth = 0.5, alpha = 0.5) +
  geom_line(data = ~subset(.x, class %in% c("simulated sample paths", "true sample paths", 
                                            "estimated sample paths")), stat = "smooth", 
            span = 0.25, se = FALSE, alpha = 0.5, size = 0.5) +
  scale_x_continuous(name = "Time", breaks = seq(0, tau, by = 0.2)) +
  scale_y_continuous(name = "Ornstein-Uhlenbeck process", breaks = seq(-8, 6, by = 2)) +
  guides(color = 'none') +
  facet_wrap(vars(class)) +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
