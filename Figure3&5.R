library(ggplot2)
source("code/ldm.R")
source("code/kerFctn.R")
sapply(paste0('mcfda/', list.files('mcfda')), source)

load("data/bmd.RData")
n1 <- length(Lt1)# 153
n2 <- length(Lt2)# 127
M <- 100
bw1 <- c(0.1, 1)
bw2 <- bw1
agef <- c(10.1, 11:24)
z0f <- c(0.778, 10.1)
agem <- 9:24
z0m <- c(0.642, 9)
bmf <- matrix(rnorm(M * length(agef)), nrow = M)
fitf <- ldm(Ly1, Lt1, z0f, agef, optns = list(M = M, bw1 = bw1, bw2 = bw2, 
                                              regular = FALSE, bm = bmf))
fitspf <- sp(Ly1, Lt1, bw = 1, domain = c(10, 24), newt = agef, x0 = z0f[1], 
             M = M, bm = bmf)# Lin and Wang (2022) JASA
bmm <- matrix(rnorm(M * length(agem)), nrow = M)
fitm <- ldm(Ly2, Lt2, z0m, agem, optns = list(M = M, bw1 = bw1, bw2 = bw2, 
                                              regular = FALSE, bm = bmm))
fitspm <- sp(Ly2, Lt2, bw = 1, domain = c(9, 24), newt = agem, x0 = z0m[1], 
             M = M, bm = bmm)# Lin and Wang (2022) JASA
load("data/bmdplot.RData")
df <- data.frame(Age = c(unlist(Lt1), unlist(Lt2), rep(agef, M), 
                         rep(agem, M), rep(agef, 3), rep(agem, 3), 10.1, 9), 
                 sbmd = c(unlist(Ly1), unlist(Ly2), t(fitf$path), t(fitm$path), 
                          t(apply(fitf$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                          t(apply(fitm$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                          0.778, 0.642), 
                 id = rep(1:(n1 + n2 + M * 2 + 6 + 2), 
                          c(sapply(Ly1, length), sapply(Ly2, length), 
                            rep(length(agef), M), rep(length(agem), M), 
                            rep(length(agef), 3), rep(length(agem), 3), 1, 1)), 
                 group = factor(rep(1:2, 
                                    c(sum(sapply(Ly1, length)) + sum(sapply(Ly2, length)), 
                                      M * length(agef) + M * length(agem) + 3 * length(agef) + 
                                        3 * length(agem) + 1 + 1)), levels = c(1, 2), 
                                labels = c('observed density snippets', 'estimated density curves')), 
                 sex = rep(c('female', 'male', 'female', 'male', 
                             'female', 'male', 'female', 'male'), 
                           c(sum(sapply(Ly1, length)), sum(sapply(Ly2, length)), 
                             M * length(agef), M * length(agem), 3 * length(agef), 3 * length(agem), 1, 1)), 
                 pred = rep(1:4, c(sum(sapply(Ly1, length)) + sum(sapply(Ly2, length)), 
                                   M * length(agef) + M * length(agem), 3 * length(agef) + 
                                     3 * length(agem), 1 + 1)))
pdf(file = "latex/img/bmd.pdf", width = 12, height = 8)
ggplot(data = df,
       aes(x = Age, y = sbmd, group = id, color = sex)
) +
  geom_line(data = ~subset(.x, pred == 1)) +
  geom_line(data = ~subset(.x, pred == 2), stat = "smooth", span = 0.4, 
            se = FALSE, alpha = 0.3) +
  geom_line(data = ~subset(.x, pred == 3), stat = "smooth", span = 0.4, 
            se = FALSE, color = 'black', linewidth = 0.7, linetype = "longdash") +
  geom_point(data = ~subset(.x, pred == 4), size = 2, color = 'black') +
  geom_text(data = ~subset(.x, pred == 4), aes(x = Age + 0.5, y = sbmd - 0.05, 
                                               label = sbmd), color = 'black') +
  facet_grid(rows = vars(sex), cols = vars(group), scales = 'free_x') +
  scale_color_discrete(guide = 'none') +
  scale_x_continuous(name = 'Age (year)', breaks = seq(8, 26, 2)) +
  scale_y_continuous(name = 'Spinal bone mineral density', 
                     breaks = seq(0.5, 1.5, 0.25), limits = c(0.5, 1.5)) + 
  theme_bw() +
  theme(text = element_text(size=20))
dev.off()

nodef <- read.table('node/bmd/bmdf.txt')
nodem <- read.table('node/bmd/bmdm.txt')
nodef <- nodef[sample(1:nrow(nodef), M), ]
nodem <- nodem[sample(1:nrow(nodem), M), ]
dfc <- data.frame(Age = c(rep(agef, M), rep(agem, M), rep(agef, 3), rep(agem, 3), 10.1, 9, 
                          rep(agef, M), rep(agem, M), rep(agef, 3), rep(agem, 3), 10.1, 9,
                          rep(agef, M), rep(agem, M), rep(agef, 3), rep(agem, 3), 10.1, 9), 
                  sbmd = c(t(fitspf$path), t(fitspm$path), 
                           t(apply(fitspf$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                           t(apply(fitspm$path, 2, quantile, c(0.05, 0.5, 0.95))), 0.778, 0.642,
                           unlist(t(nodef)), unlist(t(nodem)),
                           t(apply(nodef, 2, quantile, c(0.05, 0.5, 0.95))), 
                           t(apply(nodem, 2, quantile, c(0.05, 0.5, 0.95))), 0.778, 0.642,
                           t(fitf$path), t(fitm$path), 
                           t(apply(fitf$path, 2, quantile, c(0.05, 0.5, 0.95))), 
                           t(apply(fitm$path, 2, quantile, c(0.05, 0.5, 0.95))), 0.778, 0.642), 
                  id = rep(1:(M * 6 + 18 + 6), 
                           c(rep(length(agef), M), rep(length(agem), M), 
                             rep(length(agef), 3), rep(length(agem), 3), 1, 1, 
                             rep(length(agef), M), rep(length(agem), M), 
                             rep(length(agef), 3), rep(length(agem), 3), 1, 1,
                             rep(length(agef), M), rep(length(agem), M), 
                             rep(length(agef), 3), rep(length(agem), 3), 1, 1)), 
                  group = factor(rep(c('LW', 'Neural ODE', 'DM'), each = M * length(agef) + 
                                       M * length(agem) + 3 * length(agef) + 
                                       3 * length(agem) + 1 + 1), levels = c('DM', 'LW', 'Neural ODE')), 
                  sex = rep(c('female', 'male', 'female', 'male', 'female', 'male', 
                              'female', 'male', 'female', 'male', 'female', 'male', 
                              'female', 'male', 'female', 'male', 'female', 'male'), 
                            c(M * length(agef), M * length(agem), 3 * length(agef), 
                              3 * length(agem), 1, 1, M * length(agef), M * length(agem), 
                              3 * length(agef), 3 * length(agem), 1, 1, M * length(agef), 
                              M * length(agem), 3 * length(agef), 3 * length(agem), 1, 1)), 
                  pred = rep(rep(1:3, c(M * length(agef) + M * length(agem), 
                                        3 * length(agef) + 3 * length(agem), 1 + 1)), 3))
pdf(file = "latex/img/bmdc.pdf", width = 13.5, height = 6)
ggplot(data = dfc,
       aes(x = Age, y = sbmd, group = id, color = sex)
) +
  geom_line(data = ~subset(.x, pred == 1), stat = "smooth", span = 0.4, 
            se = FALSE, alpha = 0.3) +
  geom_line(data = ~subset(.x, pred == 2), stat = "smooth", span = 0.4, 
            se = FALSE, color = 'black', linewidth = 0.7, linetype = "longdash") +
  geom_point(data = ~subset(.x, pred == 3), size = 2, color = 'black') +
  geom_text(data = ~subset(.x, pred == 3), aes(x = Age + 0.5, y = sbmd - 0.05, 
                                               label = sbmd), color = 'black') +
  facet_grid(rows = vars(sex), cols = vars(group), scales = 'free_x') +
  scale_color_discrete(guide = 'none') +
  scale_x_continuous(name = 'Age (year)', breaks = seq(8, 26, 2)) +
  scale_y_continuous(name = 'Spinal bone mineral density', 
                     breaks = seq(0.5, 1.5, 0.25)) +
  theme_bw() +
  theme(text = element_text(size=20))
dev.off()
