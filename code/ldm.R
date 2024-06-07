#' @title Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets
#' @description Dynamic modeling of sparse longitudinal data and functional
#' snippets with stochastic differential equations.
#' @param Ly a list of \eqn{n} vectors containing the observed values
#' for each individual.
#' @param Lt a list of \eqn{n} vectors containing the observation time points
#' for each individual corresponding to \code{Ly}. Each vector should be
#' sorted in ascending order.
#' @param z0 initial condition (starting value and starting time) for
#' the underlying stochastic process.
#' @param tp a vector of length \eqn{K} of discretized time points at which
#' the stochastic process is evaluated
#' @param optns a list of options control parameters specified by
#' \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{M}{a scalar holding the number of Monte Carlo simulations to run. Default is 100.}
#' \item{regular}{whether to assume regular observation time points
#' (time spacing is the same for all individuals). Default is TRUE.}
#' \item{kernel}{smoothing kernel choice, common for mean and variance.
#' Available options are 'gauss', 'rect', 'epan', 'gausvar' and 'quar'.
#' Default is 'gauss'.}
#' \item{bw1}{bandwidth for conditional mean estimation, if not entered
#' it would be chosen from cross validation.}
#' \item{bw2}{bandwidth for conditional variance estimation, if not entered
#' it would be chosen from cross validation.}
#' \item{bm}{Brownian motion path used.}
#' @return A \code{dm} object --- a list containing the following fields:
#' \item{path}{a \eqn{M} by \eqn{K} matrix holding the \eqn{M} estimated sample
#' paths of the stochastic process at the \eqn{K} discretized time points.}
#' \item{Ly}{the original \code{Ly} used.}
#' \item{Lt}{the original \code{Lt} used.}
#' \item{z0}{the initial condition used.}
#' \item{t}{the discretized time points used.}
#' \item{optns}{the control options used.}
#' @examples
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2023. Dynamic Modeling of Sparse
#' Longitudinal Data and Functional Snippets With Stochastic Differential
#' Equations. arXiv preprint arXiv:2306.10221.}
#' }
#' @export

ldm <- function(Ly = NULL,
                Lt = NULL,
                z0 = NULL,
                tp = NULL,
                optns = list()) {
  if (is.null(Ly) | is.null(Lt) | is.null(z0) | is.null(tp)) {
    stop("require the input of Ly, Lt, z0, and tp")
  }
  if (!is.list(Ly) | !is.list(Lt)) {
    stop("Ly and Lt must be lists")
  }
  if (length(Ly) != length(Lt)) {
    stop("Ly and Lt must have the same length")
  }
  if (any(sapply(Ly, length) - sapply(Lt, length))) {
    stop("the length of each individual's observation and time points must be the same")
  }
  if (is.null(optns$M)) {
    optns$M <- 100
  }
  if(is.vector(z0)) {
    if(z0[2] != tp[1]) {
      stop("the starting time must be the same as the first time point of the time grid")
    }
    z0 <- matrix(rep(z0, each = optns$M), ncol = 2)
  } else {
    if(any(unique(z0[, 2]) != tp[1])) {
      stop("the starting time must be the same as the first time point of the time grid")
    }
  }
  if (is.null(optns$regular)) {
    optns$regular <- TRUE
  }
  if (is.null(optns$kernel)) {
    optns$kernel <- "gauss"
  }
  K <- length(tp)# number of discretized time points
  if (is.null(optns$bm)) {
    optns$bm <- matrix(rnorm((K - 1) * optns$M), nrow = optns$M)
  } else if (is.vector(optns$bm)) {
    optns$bm <- t(optns$bm)
  }
  n <- length(Ly)# number of individuals
  y <- NULL
  z <- NULL
  for (i in 1:n) {
    Ni <- length(Ly[[i]])
    y <- c(y, Ly[[i]][-1])
    if (optns$regular) {
      z <- rbind(z, cbind(Ly[[i]][-Ni], Lt[[i]][-Ni]))
    } else {
      z <- rbind(z, cbind(Ly[[i]][-Ni], Lt[[i]][-Ni], Lt[[i]][-1]))
    }
  }
  p <- ncol(z)# number of covariates
  colnames(z) <- c('y1', 't', 's')[1:p]
  p0 <- min(p, 2)
  zt <- z[, 1:min(p, 2)]
  N <- length(y)# number of pairs
  
  # select kernel
  kern <- kerFctn(optns$kernel)
  KF <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * kern(x[i] / h[i])
    }
    return(as.numeric(k))
  }
  
  # choose bandwidth by cross-validation
  if (is.null(optns$bw1)) {
    hs <- matrix(0, p0, 20)
    for (l in 1:p0) {
      hs[l, ] <- exp(seq(
        from = log(N ^ (-1 / (1 + p0)) * (max(zt[, l]) - min(zt[, l])) / 10),
        to = log(5 * N ^ (-1 / (1 + p0)) * (max(zt[, l]) - min(zt[, l]))),
        length.out = 20
      ))
    }
    
    cv <- rep(0, 20 ^ p0)
    for (k in 0:(20 ^ p0 - 1)) {
      h <- array(0, p0)
      for (l in 1:p0) {
        kl <- floor((k %% (20 ^ l)) / (20 ^ (l - 1))) + 1
        h[l] <- hs[l, kl]
      }
      for (j in 1:N) {
        a <- zt[j, ]
        if (p0 > 1) {
          mu1 <- rowMeans(apply(as.matrix(zt[-j, ]), 1, function(zi)
            KF(zi - a, h) * (zi - a)))
          mu2 <- matrix(rowMeans(apply(as.matrix(zt[-j, ]), 1, function(zi)
            KF(zi - a, h) * ((zi - a) %*% t(zi - a)))), ncol = p0)# p0^2
        } else {
          mu1 <- mean(apply(as.matrix(zt[-j, ]), 1, function(zi)
            KF(zi - a, h) * (zi - a)))
          mu2 <- mean(apply(as.matrix(zt[-j, ]), 1, function(zi)
            KF(zi - a, h) * ((zi - a) %*% t(zi - a))))
        }
        skip <- FALSE
        tryCatch(
          solve(mu2),
          error = function(e)
            skip <<- TRUE
        )
        if (skip) {
          cv[k + 1] <- Inf
          break
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(as.matrix(zt[-j, ]), 1, function(zi) {
          KF(zi - a, h) * (1 - wc %*% (zi - a))
        })# weight
        yj <- weighted.mean(y[-j], w)
        cv[k + 1] <- cv[k + 1] + (yj - y[j]) ^ 2 / N
      }
    }
    
    bwi <- which.min(cv)
    optns$bw1 <- array(0, p0)
    for (l in 1:p0) {
      kl <- floor((bwi %% (20 ^ l)) / (20 ^ (l - 1))) + 1
      optns$bw1[l] <- hs[l, kl]
    }
  } else {
    if (sum(optns$bw1 <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw1) != p0) {
      stop("dimension of bandwidth should be 2")
    }
  }
  if (is.null(optns$bw2)) {
    optns$bw2 <- optns$bw1
  } else {
    if (sum(optns$bw2 <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw2) != p0) {
      stop("dimension of bandwidth should be 2")
    }
  }
  if(!optns$regular) {
    optns$bw1 <- optns$bw1[c(1, 2, 2)]
    optns$bw2 <- optns$bw2[c(1, 2, 2)]
  }
  
  cm <- rep(0, N)
  for (i in 1:N) {
    a <- z[i, ]
    if (p > 1) {
      mu1 <- rowMeans(apply(z, 1, function(zi)
        KF(zi - a, optns$bw1) * (zi - a)))
      mu2 <- matrix(rowMeans(apply(z, 1, function(zi)
        KF(zi - a, optns$bw1) * ((zi - a) %*% t(zi - a)))), ncol = p)# p^2
    } else {
      mu1 <- mean(apply(z, 1, function(zi)
        KF(zi - a, optns$bw1) * (zi - a)))
      mu2 <- mean(apply(z, 1, function(zi)
        KF(zi - a, optns$bw1) * ((zi - a) %*% t(zi - a))))
    }
    skip <- FALSE
    tryCatch(
      solve(mu2),
      error = function(e)
        skip <<- TRUE
    )
    if (skip) {
      warning("system is exactly singular when computing fitted conditional mean")
      cm[i] <- y[i]
      # cm[i] <- NA
      break
    }
    wc <- t(mu1) %*% solve(mu2) # 1 by p
    w <- apply(z, 1, function(zi) {
      KF(zi - a, optns$bw1) * (1 - wc %*% (zi - a))
    })# weight
    cm[i] <- weighted.mean(y, w)
  }
  rcov <- (y - cm) ^ 2
  
  path <- matrix(0, nrow = optns$M, ncol = K)
  path[, 1] <- z0[, 1]
  for (q in 1:optns$M) {
    for (i in 2:K) {
      # 1st element of t is the starting time
      a <- c(path[q, i - 1], tp[i - 1])
      if (!optns$regular) {
        a <- c(a, tp[i])
      }
      if (p > 1) {
        mu1 <- rowMeans(apply(z, 1, function(zi)
          KF(zi - a, optns$bw1) * (zi - a)))
        mu2 <- matrix(rowMeans(apply(z, 1, function(zi)
          KF(zi - a, optns$bw1) * ((zi - a) %*% t(zi - a)))), ncol = p)# p^2
      } else {
        mu1 <- mean(apply(z, 1, function(zi)
          KF(zi - a, optns$bw1) * (zi - a)))
        mu2 <- mean(apply(z, 1, function(zi)
          KF(zi - a, optns$bw1) * ((zi - a) %*% t(zi - a))))
      }
      skip <- FALSE
      tryCatch(
        solve(mu2),
        error = function(e)
          skip <<- TRUE
      )
      if (skip) {
        warning("system is exactly singular when generating simulated trajectory")
        path[q, i:K] <- rep(NA, K - i + 1)
        break
      }
      wc <- t(mu1) %*% solve(mu2)
      w <- apply(z, 1, function(zi) {
        KF(zi - a, optns$bw1) * (1 - wc %*% (zi - a))
      })# weight
      cmi <- weighted.mean(y, w)# conditional mean
      cci <- max(weighted.mean(rcov, w), 0)# conditional variance
      path[q, i] <- optns$bm[q, i - 1] * sqrt(cci) + cmi
    }
  }
  res <- list(
    path = path,
    Ly = Ly,
    Lt = Lt,
    z0 = z0,
    tp = tp,
    optns = optns
  )
  class(res) <- "dm"
  res
}