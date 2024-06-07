#' Estimate Mean Function
#' @param t a list of vectors (for irregular design) or a vector (for regular design) containing time points of observations for each individual. Each vector should be in ascending order
#' @param y a list of vectors (for irregular design) or a matrix (for regular design) containing the observed values at \code{t}. If it is a matrix, the columns correspond to the time points in the vector \code{t}
#' @param newt  a list of vectors or a vector containing time points of observations to be evaluated. If NULL, then newt is treated as t
#' @param method estimation method, 'PACE' or 'FOURIER'
#' @param tuning tuning method to select possible tuning parameters
#' @param ... other parameters required depending on the \code{method} and \code{tuning}; see details
#' @details
#'     \itemize{
#'         \item{When \code{method='PACE'}, additional parameters are}
#'         \describe{
#'             \item{\code{kernel}}{kernel type; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic", "sigmoid" and "silverman"; see https://en.wikipedia.org/wiki/Kernel_(statistics) for more details.}
#'             \item{\code{deg}}{degree of the local polynomial regression; currently only \code{deg=1} is supported.}
#'             \item{\code{bw}}{bandwidth}
#'         }
#'         \item{When \code{method='FOURIER'}, additional parameters are}
#'         \describe{
#'             \item{\code{q}}{number of basis functions; if \code{NULL} then selected by \code{tuning} method}
#'             \item{\code{rho}}{roughness penalty parameter; if \code{NULL} then selected by \code{tuning} method}
#'             \item{\code{ext}}{extension margin of Fourier extension; if \code{NULL} then selected by \code{tuning} method}
#'             \item{\code{domain}}{time domain; if \code{NULL} then estimated by \code{(min(t),max(t))}}
#'         }
#'     }
#'
#' @return an object of the class 'meanfunc' containing necessary information to predict the mean function
#' \itemize{
#'     \item{When \code{method='PACE'}, additional parameters are}
#'         \describe{
#'             \item{\code{fitted}}{fitted value at \code{newt}}
#'             \item{\code{bw}}{selected bandwidth by \code{tuning} method if \code{NULL} is the input for \code{bw}.}
#'         }
#'     \item{When \code{method='FOURIER'}, additional parameters are}
#'         \describe{
#'             \item{\code{fitted}}{fitted value at \code{newt}}
#'             \item{\code{q}}{selected \code{q} if \code{NULL} is the input}
#'             \item{\code{rho}}{selected \code{rho} if \code{NULL} is the input}
#'             \item{\code{ext}}{selected \code{ext} if \code{NULL} is the input}
#'             \item{\code{bhat}}{estimated coefficients}
#'         }
#' }
#' 
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lin2020}{mcfda}
#' 
#' \insertRef{Yao2005}{mcfda}
#' 
#' @examples
#' mu <- function(s) sin(2*pi*s)
#' D <- synfd::sparse.fd(mu=mu, X=synfd::gaussian.process(), n=100, m=5)
#' mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='PACE',
#'                 tuning='cv',weig=NULL,kernel='gauss',deg=1)
#' # equivalent to
#' # mu.obj <- meanfunc(D$t,D$y)
#'
#' # plot the object
#' plot(mu.obj)
#' @export meanfunc

meanfunc <- function(Ly,Lt,bw=NULL,domain=NULL,newt=NULL,kernel='gauss')
{
    if(is.null(domain))
    {
        domain <- c(0,0)
        domain[1] <- min(unlist(Lt))
        domain[2] <- max(unlist(Lt))
        domain <- c(domain[1]-0.01*(domain[2]-domain[1]),
                    domain[2]+0.01*(domain[2]-domain[1]))
    }
    
    n <- length(Lt)
    x <- unlist(Lt)
    y <- unlist(Ly)
    ord <- sort(x,index.return=T)$ix
    x <- x[ord]
    y <- y[ord]
    
    if(is.null(bw)) bw <- bw.lp1D(x,y,kernel=kernel,method='cv',K=5,H=NULL)

    R <- list(bw=bw,x=x,y=y,n=n,domain=domain,kernel=kernel,yend=c(NULL,NULL))
    L0 <- domain[2]-domain[1]
    yend <- predictmf(R,c(domain[1]+L0/100,domain[2]-L0/100))
    R$yend <- yend
    if(!is.null(newt)) {
      R$fitted <- predictmf(R,newt)
    }
    return(R)
}


#' predict mean functions at new locations
#' @param meanfunc.obj the object obtained by calling \code{mean.func}
#' @param newt a vector or a list of vectors of real numbers
#' @return the estimated mean function evaluated at \code{newt}. It has the same format of \code{newt}
#' @export
predictmf <- function(meanfunc.obj,newt)
{
    pred <- function(newt) # newt must be a vector
    {
        idxl <- newt < meanfunc.obj$domain[1]
        idxu <- newt > meanfunc.obj$domain[2]
        idx <- (!idxl) & (!idxu)
        
        newt0 <- newt[idx]
        ord <- sort(newt0,index.return=T)$ix
        
        tmp <- rep(Inf,length(newt0))
        
        tmp[ord] <- fdapace::Lwls1D(bw = meanfunc.obj$bw, 
                                    kernel_type = meanfunc.obj$kernel,
                                    xin = meanfunc.obj$x,
                                    yin = meanfunc.obj$y,
                                    xout = newt0[ord])
        

        # tmp <- loclin1D(x=meanfunc.obj$x,
        #                        y=meanfunc.obj$y,
        #                        newx=newt[idx],
        #                        bw=meanfunc.obj$bw,
        #                        weig=meanfunc.obj$weig,
        #                        kernel=meanfunc.obj$kernel,
        #                        deg=meanfunc.obj$deg)$fitted

        yhat <- rep(0,length(newt))
        yhat[idx] <- tmp
        yhat[idxl] <- meanfunc.obj$yend[1]
        yhat[idxu] <- meanfunc.obj$yend[2]
        return(yhat)
    }

    if(is.list(newt))
    {
        mi <- lapply(newt,length)
        newt <- unlist(newt)
        fitted <- pred(newt)

        cm <- c(0,cumsum(mi))

        R <- sapply(1:length(mi),function(i){
            res <- list()
            res[[1]] <- fitted[(cm[i]+1):cm[i+1]]
            res
        })
        return(R)
    }
    else if(is.vector(newt))
    {
        return(pred(newt))
    }
    else stop('newt must be a vector or a list of vectors of real numbers')
}