#' @title Create covariance function for a fitted OU or OUF ctmm model object
#' @description A function is created to evaluate the covariance function of the fitted
#' OUF movement model
#' @param x A `ctmm` object created by a call to \code{\link[ctmm]{ctmm.fit}} or \code{\link[ctmm]{ctmm.select}}
#' @details The function returns a function to evaluate the covariance of the fitted
#' Ornstein-Ulenbeck Foraging movement model. The returned function has 2 arguments: (1) \code{s}
#' and (2) \code{t} both vectors of times to evaluate the covariance function of the fitted
#' OUF model.
#'@export
#'@author Devin S. Johnson
#'
ouf_corfun <- function(x){
  sss <- summary(x)
  is_iid <- sss$name %in%  c("IID anisotropic", "IID")
  is_ou <- sss$name %in%  c("OU anisotropic", "OU")
  is_ouf <- sss$name %in%  c("OUF anisotropic", "OUF", "OUf anisotropic", "OUf")
  if(!(is_ou | is_ouf | is_iid)) stop("'",paste0(sss$name), "' model not currently supported with this function.")
  if(is_iid){
    foo <- function(s,t){1.0*(s==t)}
    attr(foo, "tau") <- c(0,0)
  } else if(is_ou){
    tau_p <- x$tau[[1]]
    if(tau_p<1.0e-8){
      warning("OU model seems to be approx. IID, using IID correlation function!")
      foo <- function(s,t){1.0*(s==t)}
      attr(foo, "tau") <- c(0,0)
    } else{
      foo <- function(s, t){exp(-abs(t-s)/tau_p)}
    }
    attr(foo, "tau") <- c(tau_p,0)
  } else{
    tau_p <- x$tau[[1]]
    tau_v <- x$tau[[2]]
    if(tau_v<1.0e-8 & tau_p>1.0e-8){
      warning("OUF(f) model seems to be approx. OU, using OU correlation function!")
      foo <- function(s, t){exp(-abs(t-s)/tau_p)}
      attr(foo, "tau") <- c(tau_p, 0)
    } else if(tau_v<1.0e-8 & tau_p<1.0e-8){
      warning("OUF(f) model seems to be approx. IID, using IID correlation function!")
      foo <- function(s,t){1.0*(s==t)}
      attr(foo, "tau") <- c(0,0)
    } else{
      stop("OUF seems to be occilitory, i.e, tau_v > tau_p.")
    }
  }
  return(foo)
}



#' @title Calculate correlation matrix for a set of times from a `ctmm` covariance function
#' @description Using a correlation function created by \code{\link{ouf_corfun}}
#' from a fitted OUF related model a covariance (correlation) matrix is created for observations
#' at the user provided times.
#' @param x A \code{ctmm} model object from a call to \code{\link[ctmm]{ctmm.fit}}
#' or `\link[ctmm]{ctmm.select}`.
#' @param times A vector of POSIX times at which the covariance matrix will be constructed.
#' A `telemetry` data object will also work.
#' @param inverse Logical. Should the inverse covariance matrix be returned. Defaults to `FALSE`
#' @author Devin S. Johnson
#' @export
ouf_covmat <- function(x, times, inverse=FALSE){
  if(!inherits(x, "ctmm")) stop("'x' must be a 'ctmm' object.")
  if(inherits(times, "telemetry")) times <- times$timestamp
  if(!inherits(times, "POSIXct")) stop("'times' must be a POSIX vector or 'telemetry' data object.")
  times <- as.numeric(times)
  times <- times-min(times)
  cf <- ouf_corfun(x)
  Sig0 <- x$sigma
  R <- outer(times, times, cf)
  if(inverse){
    out <- kronecker(solve(R), solve(Sig0))
  } else{
    out <- kronecker(R, Sig0)
  }
  return(out)
}
