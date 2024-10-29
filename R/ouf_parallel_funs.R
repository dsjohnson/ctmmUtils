#' @title Parallel fitting of ctmm OUF family models
#' @param tdata A list of telemetry data objects.
#' @param ... Additional arguments (besides `data` and `CTMM`!) passed to `\link[ctmm]{ctmm.select}`
#' @importFrom progressr progressor
#' @importFrom ctmm ctmm.guess ctmm.select ctmm
#' @import foreach doFuture
#' @export
ctmm_select_parallel <- function(tdata, ...){
  i <- NULL
  xargs <- list(...)
  arg_verbose <- ifelse(!is.null(xargs$verbose), xargs$verbose, FALSE)
  if(arg_verbose & add_ess){
    add_ess <- FALSE
    warning("Currently cannot calculate ESS for multiple models. Skipping ESS calculation.")
  }
  progressr::handlers(global = TRUE)
  p <- progressr::progressor(length(tdata))
  out <- foreach(i=1:length(tdata),
                 .options.future = list(seed = TRUE),.errorhandling = "pass") %dofuture% {
                   guess <- ctmm.guess(tdata[[i]], ctmm(error=TRUE), interactive=FALSE)
                   fit <- ctmm.select(tdata[[i]], guess,...)
                   # if(add_ess){
                   #   sum_fit <- summary(fit)
                   #   include <- c("IID","OU","OUF","IID anisotropic","OU anisotropic","OUF anisotropic")
                   #   fit <- fit[c(1:nrow(sum_fit))[rownames(sum_fit) %in% include]][[1]]
                   #   ess <- ouf_ess(fit, tdata[[i]])
                   #   fit <- list(fit=fit, ess=ess)
                   # }
                   p()
                   fit
                 }
  if(!add_ess){
    chk <- all(sapply(out, class)=="ctmm")
    if(!chk) warning("There appears to be some fitting problems check the list of fitted models!")
  } else{
    chk <- all(sapply(out, \(x) class(x$fit))=="ctmm")
    if(!chk) warning("There appears to be some model fitting problems!")
    chk2 <- all(sapply(out, \(x) inherits(x$ess,"oufESS")))
    if(!chk2) warning("There appears to be some ESS calculation problems!")
  }

  return(out)
}

#' @title Parallel ESS calculation of ctmm OUF family models
#' @param fits A list of fitted `ctmm` OUF family models
#' @param tdata A list of telemetry data objects.
#' @param ... Additional arguments (besides `data` and `CTMM`!) passed to `\link[ctmm]{ctmm.select}`
#' @importFrom progressr progressor
#' @import foreach doFuture
#' @export
ouf_ess_parallel <- function(fits, tdata, ...){
  i <- NULL
  progressr::handlers(global = TRUE)
  p <- progressr::progressor(length(fits))
  out <- foreach(i=1:length(fits),
                 .options.future = list(seed = TRUE),
                 .errorhandling = "pass") %dofuture% {
                   ess <- ctmmUtils::ouf_ess(fits[[i]], tdata[[i]])
                   p()
                   ess
                 }
  chk <- all(sapply(ess, inherits, "oufESS"))
  if(!chk) warning("There appears to be some issues, check the list of ess calculations!")
  return(out)
}

