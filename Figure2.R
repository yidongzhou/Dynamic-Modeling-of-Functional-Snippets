library(ggplot2)
source("code/ldm.R")
source("code/kerFctn.R")

load("data/ngd.RData")
n1 <- length(Lt1)# 87
n2 <- length(Lt2)# 96
M <- 100
bw1 <- c(5, 10)
bw2 <- bw1
agef <- seq(4, 72, 4)
z0f <- c(52.9, 4)
fitf <- ldm(Ly = Ly1, Lt = Lt1, z0 = z0f, tp = agef, optns = list(M = M, bw1 = bw1, bw2 = bw2))
bw1 <- c(8, 2.5)
bw2 <- bw1
agem <- seq(12, 72, 4)
z0m <- c(63.0, 12)# 20: 65.1; 36: 80.3
fitm <- ldm(Ly = Ly2, Lt = Lt2, z0 = z0m, tp = agem, optns = list(M = M, bw1 = bw1, bw2 = bw2))
# load("data/ngdplot.RData")
df <- data.frame(Age = c(unlist(Lt1), unlist(Lt2), rep(agef, M), rep(agem, M), 
                         rep(agef, 3), rep(agem, 3), 4, 12, 20), 
                 Height = c(unlist(Ly1), unlist(Ly2), t(fitf$path), t(fitm$path), 
                            t(apply(fitf$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                            t(apply(fitm$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                            52.9, 63, 65.1), 
                 id = rep(1:(n1 + n2 + M * 2 + 6 + 2), 
                          c(sapply(Ly1, length), sapply(Ly2, length), 
                            rep(length(agef), M), rep(length(agem), M), 
                            rep(length(agef), 3), rep(length(agem), 3), 1, 2)), 
                 group = factor(rep(1:2, 
                                    c(sum(sapply(Ly1, length)) + sum(sapply(Ly2, length)), 
                                      M * length(agef) + M * length(agem) + 3 * length(agef) + 
                                        3 * length(agem) + 1 + 2)), levels = c(1, 2), 
                                labels = c('observed growth snippets', 'estimated growth curves')), 
                 sex = rep(c('female', 'male', 'female', 'male', 
                             'female', 'male', 'female', 'male'), 
                           c(sum(sapply(Ly1, length)), sum(sapply(Ly2, length)), 
                             M * length(agef), M * length(agem), 
                             3 * length(agef), 3 * length(agem), 1, 2)), 
                 pred = rep(1:4, c(sum(sapply(Ly1, length)) + 
                                     sum(sapply(Ly2, length)), 
                                   M * length(agef) + M * length(agem), 
                                   3 * length(agef) + 3 * length(agem), 1 + 2)))
pdf(file = "latex/img/ngd.pdf", width = 12, height = 8)
ggplot(data = df,
       aes(x = Age, y = Height, group = id, color = sex)
) +
  geom_line(data = ~subset(.x, pred == 1)) +
  geom_line(data = ~subset(.x, pred == 2), stat = "smooth", span = 0.4, 
            se = FALSE, alpha = 0.3) +
  geom_line(data = ~subset(.x, pred == 3), stat = "smooth", span = 0.4, 
            se = FALSE, color = 'black', linewidth = 0.7, linetype = "longdash") +
  geom_point(data = ~subset(.x, pred == 4), size = 2, color = 'black') +
  geom_text(data = ~subset(.x, pred == 4), aes(x = Age + 2, y = Height - 2, 
                                               label = Height), color = 'black') +
  facet_grid(rows = vars(sex), cols = vars(group), scales = 'free_x') +
  scale_color_discrete(guide = 'none') +
  # scale_alpha_continuous(guide = 'none') +
  scale_x_continuous(name = 'Age (month)', breaks = seq(0, 72, 8)) +
  scale_y_continuous(name = 'Height (cm)', breaks = seq(50, 120, 10), 
                     limits = c(50, 120)) + 
  theme_bw() +
  theme(text = element_text(size=20))
dev.off()
