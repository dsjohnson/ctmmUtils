#' @title Parallel fitting of ctmm OUF family models
#' @param tdata A list of telemetry data objects.
#' @param ... Additional arguments (besides `data` and `CTMM`!) passed to `\link[ctmm]{ctmm.select}`
#' @importFrom progressr progressor
#' @importFrom ctmm ctmm.guess ctmm.select
#' @import foreach doFuture

#' @export
ctmm_select_parallel <- function(tdata, ...){
  i <- NULL
  progressr::handlers(global = TRUE)
  p <- progressr::progressor(length(tdata))
  out <- foreach(i=1:length(tdata),
                 .options.future = list(seed = TRUE),
                 .errorhandling = "pass") %dofuture% {
                   guess <- ctmm.guess(tdata[[i]], interactive=FALSE)
                   suppressWarnings(fit <- ctmm.select(tdata[[i]], guess, ...))
                   p()
                   fit
                 }
  chk <- all(sapply(out, class)=="ctmm")
  if(!chk) warning("There appears to be some fitting problems check the list of fitted models!")
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
  if(!chk) warning("There appears to be some problems, check the list of ess calculations!")
  return(out)
}

