#' @title Calculate Effective Sample Size for a Set of OUF locations
#' @description Estimates the number of independent locations in a `ctmm` data set
#' using the mutual information method of Bartoszek (2016).
#' @param x A \code{ctmm} object (See \code{\link[ctmm]{ctmm.fit}}).
#' @param times A vector of POSIX times. A `telemtry` data object will also work.
#' @details This function uses the "mutual information" effective sample size of Bartoszek (2016) to calculate the
#' equivalent number of independent animal locations. It also calculates individual contributions of each location using the
#' regression effective sample size in Bartoszek (2016). The output is a named list with `Ne` equal to the overall sample
#' and `w` is a vector of weights that sum to 1 overall. If you want the ESS value of each observation `Ne * w` will provide it.
#' @references Bartoszek, K. (2016). Phylogenetic effective sample size.
#' Journal of Theoretical Biology. 407:371-386. (See https://arxiv.org/pdf/1507.07113.pdf).
#' @author Devin S. Johnson
#' @export
#'
ouf_ess <- function(x, times){
  if(inherits(x, "list") & "ess"%in%names(x)) stop("ESS has been calculated for this model fit already. See 'x$ess'")
  if(!inherits(x, "ctmm")) stop("'x' must be a 'ctmm' object.")
  if(inherits(times, "telemetry")) times <- times$timestamp
  if(!inherits(times, "POSIXct")) stop("'times' must be a POSIX vector or 'telemetry' data object.")
  times <- as.numeric(times)
  times <- times-min(times)
  n <- length(times)
  cf <- ouf_corfun(x)
  tau <- attr(cf, "tau")
  if(tau[[1]]==0){
    out <- list(N_mi = n, N_r=n, w=rep(1/n,n))
  } else{
    Sig0 <- x$sigma
    R <- outer(times, times, cf)
    log_det_R <- as.numeric(determinant(R,logarithm = TRUE)$modulus)
    log_det_Sig0 <- as.numeric(determinant(Sig0,logarithm = TRUE)$modulus)
    n_mi <-  1 + (n-1)/log(exp(1) - log_det_R)
    Vi <- kronecker(solve(R), solve(Sig0))
    Vid <- fatdiag(Vi,2)
    w <- 1/n + (n-1)*exp(sapply(Vid, \(x) -as.numeric(determinant(x, logarithm = TRUE)$modulus)) - log_det_Sig0)/n
    n_r <- sum(w)
    w <- w/n_r
    out <- list(N_mi = n_mi, N_r=n_r, w=w)
  }
  class(out) <- c("list","oufESS")
  return(out)
}


fatdiag <- function(M,size){
  if((nrow(M)%%size != 0) |  (ncol(M)%%size != 0)) stop("nrow(M) and/or ncol(M) is not an even multiple of 'size'.")
  idx <- cbind( seq(1,nrow(M),size), seq(size,ncol(M),size))
  idx <- idx[idx[,2]>=idx[,1],]
  out <- apply(idx, 1, \(i) M[seq(i[1],i[2],1), seq(i[1],i[2],1)], simplify = FALSE)
  return(out)
  }

