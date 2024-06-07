#' estimate the variance of noise
#' @param t a list of vectors containing time points for each individual. Each vector is sorted in ascending order.
#' @param y a list of vectors containing observed values at \code{t}
#' @param h the bandwidth
#' @return the estimated variance of the measurement error
#' @export sigma2
sigma2 <- function(Ly,Lt,h=NULL,mu.hat)
{
  n <- length(Lt)

    t.min <- min(unlist(Lt))
    t.max <- max(unlist(Lt))
    
    if(is.null(h))
    {
      h <- select.sig2.bw(Ly,Lt,mu.hat)
    }
    else if(2*h >= t.max-t.min) stop('h is too large')

    AB <- sapply(1:n,function(i){
        tobs <- Lt[[i]]
        y <- Ly[[i]]
        m <- length(tobs)
        v1 <- 0
        v2 <- 0

        if(m < 2) return(c(v1,v2))
        for(j in 1:m)
        {
            for(k in 1:m)
            {
                if( (k!=j) && abs(tobs[j]-tobs[k])<h )
                {

                    v1 <- v1 + (y[j]-y[k])^2 / 2
                    v2 <- v2 + 1
                }
            }
        }
        if(v2 == 0) return(0)
        else return(v1/v2)
    })
    AB <- unlist(AB)
    nz <- sum(AB!=0)
    if(nz==0) sig2 <- 0
    else sig2 <- sum(AB)/nz

    return(sig2)
}

select.sig2.bw <- function(Ly,Lt,mu.hat)
{
  n <- length(Lt)
  
  t.min <- min(unlist(Lt))
  t.max <- max(unlist(Lt))
  
  delta <- max(sapply(Lt,function(ts){max(ts)-min(ts)}))
  m <- mean(sapply(Lt,function(ts){length(ts)}))
  M <- n * m^2
  #h <- delta * (M^(-1/5)) / 7
  
  tmp <- lapply(1:length(Lt),function(i){
    rr <- Ly[[i]] - mu.hat[[i]]
    rr^2
  })
  vn <- sqrt(mean(unlist(tmp)))
  
  h <- 0.29 * delta * vn * (M^(-1/5))
  
  max.it <- 1000
  it <- 0
  while(it < max.it) # h0 two small
  {
    it <- it + 1
    
    cnt <- sapply(1:n,function(i){
      tobs <- Lt[[i]]
      y <- Ly[[i]]
      m <- length(tobs)
      v1 <- 0
      
      if(m < 2) return(0)
      for(j in 1:m)
      {
        for(k in 1:m)
        {
          if( (k!=j) && abs(tobs[j]-tobs[k])<h )
          {
            v1 <- v1 + 1
          }
        }
      }
      return(v1)
    })
    
    cnt <- sum(cnt)
    if(cnt >= min(50,0.1*n*m*(m-1))) break
    else
    {
      h <- h*1.01
    }
  }
  
  return(h)
}