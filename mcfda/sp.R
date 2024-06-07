sp <- function(Ly, Lt, bw = NULL, domain = NULL, newt = NULL, x0, M = 100, bm = NULL) {
  mf <- meanfunc(Ly = Ly, Lt = Lt, bw = bw, domain = domain, newt = newt)
  cf <- covfunc(Ly = Ly, Lt = Lt, bw = bw, domain = domain, newt = newt, mu = mf)
  mf <- mf$fitted
  cf <- cf$fitted
  K <- length(newt)
  if (is.null(bm)) {
    bm <- matrix(rnorm((K - 1) * M), nrow = M)
  } else if (is.vector(bm)) {
    bm <- t(bm)
  }
  path <- matrix(0, nrow = M, ncol = K)
  path[, 1] <- rep(x0, M)
  for(m in 1:M) {
    for(k in 2:K) {
      if(cf[k-1, k-1] <= 0) {
        cm <- mf[k]
        cc <- cf[k, k]
      } else {
        cm <- mf[k] + cf[k, k-1] * (path[m, k-1] - mf[k-1]) / cf[k-1, k-1]
        cc <- max(cf[k, k] - cf[k, k-1] * cf[k-1, k] / cf[k-1, k-1], 0)
      }
      path[m, k] <- bm[m, k - 1] * sqrt(cc) + cm
    }
  }
  res <- list(mf = mf, sf = sqrt(diag(cf)), path = path)
  res
}