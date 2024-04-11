#' @title Parallel fitting of ctmm OUF family models
#' @param tdata A list of telemetry data objects.
#' @param add_ess Calculated "mutual information" and "regression" effective sample sizes and weights from
#' Bartoszek (2016).
#' @param ... Additional arguments (besides `data` and `CTMM`!) passed to `\link[ctmm]{ctmm.select}`
#' @importFrom progressr progressor
#' @importFrom ctmm ctmm.guess ctmm.select
#' @import foreach doFuture
#' @references Bartoszek, K. (2016). Phylogenetic effective sample size.
#' Journal of Theoretical Biology, 407, 371-386.  (See https://arxiv.org/pdf/1507.07113.pdf).

#' @export
ctmm_select_parallel <- function(tdata, add_ess=FALSE, ...){
  i <- NULL
  xargs <- list(...)
  verbose <- ifelse(!is.null(xargs$verbose), xargs$verbose, FALSE)
  if(verbose & add_ess) stop("Currently cannot set 'verbose=TRUE' and 'add_ess=TRUE'")
  progressr::handlers(global = TRUE)
  p <- progressr::progressor(length(tdata))
  out <- foreach(i=1:length(tdata),
                 .options.future = list(seed = TRUE),
                 .errorhandling = "pass") %dofuture% {
                   guess <- ctmm.guess(tdata[[i]], interactive=FALSE)
                   suppressWarnings(fit <- ctmm.select(tdata[[i]], guess, ...))
                   if(add_ess){
                     ess <- ouf_ess(fit, tdata[[i]])
                     fit <- list(fit=fit, ess=ess)
                   }
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

