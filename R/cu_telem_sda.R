#' @title Filter track for speed, distance and angle.
#' @description Applies the Freitas SDA (speed, distance, angle) filter to assess telemetry data outliers. The
#' `\link[trip]{sda}` function from the `trip` package is used for computation.
#' @param data `telemetry` data object from the `ctmm` package.
#' @param ... Additional arguments for the `link[trip]{sda}` function. The argument `smax` (in km/hr) is
#' necessary. Please see the `link[trip]{sda}` documentation.
#' @author Devin S. Johnson
#' @importFrom trip trip sda
#' @importFrom methods as is
#' @export
#'
#' @references Freitas, C., Lydersen, C., Fedak, M. A. and Kovacs, K. M. (2008),
#' A simple new algorithm to filter marine mammal Argos locations.
#' Marine Mammal Science, 24: 315?V325. doi: 10.1111/j.1748-7692.2007.00180.x

telem_sda <- function(data, ...){
  is_list <- is(data, "list")
  if(is_list){
    if(!all(sapply(data, is,  "telemetry"))) stop("'data' does not appear to be a telemetry object.")
  } else{
    if(is(data,"telemetry")) stop("'data' does not appear to be a telemetry object.")
  }

  if(!is_list){
    tr <- as(data[,c('longitude','latitude','timestamp')], "data.frame")
    tr$id <- 1
    tr <- suppressWarnings(trip(tr,  c("timestamp","id")))
    data$keep <- sda(tr,...)
    return(data)
  } else{
    for(i in 1:length(data)){
      tr <- as(data[[i]][,c('longitude','latitude','timestamp')], "data.frame")
      tr$id <- 1
      tr <- suppressWarnings(trip(tr,  c("timestamp","id")))
      data[[i]]$keep <- sda(tr,...)
    }
    return(data)
  }



}

