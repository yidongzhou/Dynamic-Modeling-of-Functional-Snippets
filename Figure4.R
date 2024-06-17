library(ggplot2)
library(fdapace)
library(fda)# Berkeley growth study data
source("code/ldm.R")
source("code/kerFctn.R")

# generate synthetic growth curves
n1 <- ncol(growth$hgtf)
n2 <- ncol(growth$hgtm)
Ly1 <- lapply(1:n1, function(i) growth$hgtf[, i])
Ly2 <- lapply(1:n2, function(i) growth$hgtm[, i])
Lt1 <- rep(list(growth$age), n1)
Lt2 <- rep(list(growth$age), n2)

fit1 <- FPCA(Ly1, Lt1, optns = list(nRegGrid = 18))
n1 <- 300
set.seed(2)
xi1 <- cbind(rnorm(n1, mean = 0, sd = sqrt(fit1$lambda[1])), 
             rnorm(n1, mean = 0, sd = sqrt(fit1$lambda[2])), 
             rnorm(n1, mean = 0, sd = sqrt(fit1$lambda[3])))
female <- rep(fit1$mu, 300) + fit1$phi[, 1:3] %*% t(xi1)

fit2 <- FPCA(Ly2, Lt2, optns = list(nRegGrid = 18))
n2 <- 300
set.seed(9)
xi2 <- cbind(rnorm(n2, mean = 0, sd = sqrt(fit2$lambda[1])), 
             rnorm(n2, mean = 0, sd = sqrt(fit2$lambda[2])), 
             rnorm(n2, mean = 0, sd = sqrt(fit2$lambda[3])))
male <- rep(fit2$mu, 300) + fit2$phi[, 1:3] %*% t(xi2)

# generate snippets from Berkeley growth study data
# randomly selecting two measurements, one year apart, for each subject
age <- 1:18
Lt1 <- list()
Ly1 <- list()
set.seed(1)
for(i in 1:n1) {
  k <- sample(length(age) - 1, 1)
  Lt1[[i]] <- age[c(k, k + 1)]
  Ly1[[i]] <- female[c(k, k+1), i]
}
Lt2 <- list()
Ly2 <- list()
set.seed(1)
for(i in 1:n2) {
  k <- sample(length(age) - 1, 1)
  Lt2[[i]] <- age[c(k, k + 1)]
  Ly2[[i]] <- male[c(k, k+1), i]
}

load("data/bgd.RData")
M <- 300
bw1 <- c(20, 0.5)
bw2 <- bw1
set.seed(1)
z0f <- cbind(sample(female[1, ], M), 1)
z0m <- cbind(sample(male[1, ], M), 1)
fitf <- ldm(Ly1, Lt1, z0f, age, optns = list(M = M, bw1 = bw1, bw2 = bw2))
fitm <- ldm(Ly2, Lt2, z0m, age, optns = list(M = M, bw1 = bw1, bw2 = bw2))
load('data/bgdplot.RData')
df <- data.frame(Age = c(rep(age, n1 + n2), unlist(Lt1), unlist(Lt2), 
                         rep(age, M + M)), 
                 Height = c(female, male, unlist(Ly1), unlist(Ly2), 
                            t(fitf$path), t(fitm$path)), 
                 id = rep(1:(n1 + n2 + n1 + n2 + M + M), 
                          c(rep(length(age), n1 + n2), 
                            rep(2, n1 + n2), rep(length(age), M + M))), 
                 group = factor(rep(1:3, c((n1 + n2) * length(age), 
                                           (n1 + n2) * 2, (M + M) * length(age))), 
                                levels = c(1, 2, 3), 
                                labels = c('synthetic growth curves', 
                                           'artificial growth snippets', 
                                           'estimated growth curves')), 
                 sex = rep(c('female', 'male', 'female', 
                             'male', 'female', 'male'), 
                           c(n1 * length(age), n2 * length(age), n1 * 2, n2 * 2, 
                             M * length(age), M * length(age))))
pdf(file = "latex/img/bgd.pdf", width = 12, height = 8)
ggplot(data = df,
       aes(x = Age, y = Height, group = id, color = factor(sex))
) +
  geom_line(data = ~subset(.x, group == "artificial growth snippets")) +
  geom_line(data = ~subset(.x, group %in% c("estimated growth curves", 
                                            "synthetic growth curves")), 
            stat = "smooth", span = 0.4, se = FALSE) +
  facet_grid(rows = vars(sex), cols = vars(group)) +
  scale_color_discrete(guide = 'none') +
  scale_alpha_continuous(guide = 'none') +
  scale_x_continuous(name = 'Age (year)', breaks = seq(0, 18, 2)) +
  scale_y_continuous(name = 'Height (cm)', breaks = seq(70, 210, 10), 
                     limits = c(65, 210)) + 
  theme_bw() +
  theme(text = element_text(size=20))
dev.off()
