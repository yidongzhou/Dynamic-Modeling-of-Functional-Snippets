#' @title Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets
#' @description This function performs dynamic modeling of sparse longitudinal 
#' data and functional snippets using stochastic differential equations. 
#' The method employs multiple linear regression to estimate both the conditional 
#' mean and the conditional variance.
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
#' \item{error}{the standard deviation of measurement errors. Default is
#' \code{NULL}, where no measurement error is assumed.}
#' \item{regular}{whether to assume regular observation time points
#' (time spacing is the same for all individuals). Default is TRUE.}
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

gdm <- function(Ly = NULL,
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
  if(z0[2] != tp[1]) {
    stop("the starting time must be the same as the first time point of the time grid")
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
  if (is.null(optns$error)) {
    optns$error <- FALSE
  }
  if (is.null(optns$regular)) {
    optns$regular <- TRUE
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
  colnames(z) <- c('y1', 't', 's')[1:ncol(z)]
  df <- data.frame(y2 = y, z)
  if (optns$error) {
    require(simex)
    fitb <- lm(y2 ~ ., data = df, x = TRUE)
    fitb <- simex(
      model = fitb,
      SIMEXvariable = 'y1',
      measurement.error = optns$error,
      asymptotic = FALSE
    )
    df$y2 <- fitb$residuals ^ 2
    fitsigma <- lm(y2 ~ ., data = df, x = TRUE)
    fitsigma <- simex(
      model = fitsigma,
      SIMEXvariable = 'y1',
      measurement.error = optns$error,
      asymptotic = FALSE
    )
  } else {
    fitb <- lm(y2 ~ ., data = df, x = TRUE)
    df$y2 <- fitb$residuals ^ 2
    fitsigma <- lm(y2 ~ ., data = df, x = TRUE)
  }
  
  path <- matrix(0, nrow = optns$M, ncol = K)
  path[, 1] <- z0[, 1]
  for (q in 1:optns$M) {
    for (i in 2:K) {
      # 1st element of tp is the starting time
      newdata <- data.frame(y1 = path[q, i - 1], t = tp[i - 1])
      if (!optns$regular) {
        newdata$s <- tp[i]
      }
      cmi <- predict(fitb, newdata)
      cci <- max(predict(fitsigma, newdata), 0)
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
