#' @title Filter track for speed, distance and angle.
#' @description Applies the Freitas SDA (speed, distance, angle) filter to assess telemetry data outliers. The
#' `\link[trip]{sda}` function from the `trip` package is used for computation.
#' @param data `telemetry` data object from the `ctmm` package.
#' @param ... Additional arguments for the `\link[trip]{sda}` function. The argument `smax` (in km/hr) is
#' necessary. Please see the `\link[trip]{sda}` documentation.
#' @param filter Logical. Remove locations determined to be outliers based on the SDA filter. If `filter=FALSE` then a
#' column named `keep` will indicate which observations should be kept, but all will remain in the data. The default is
#' `filter=TRUE`
#' @author Devin S. Johnson
#' @importFrom trip trip sda
#' @importFrom methods as is
#' @export
#'
#' @references Freitas, C., Lydersen, C., Fedak, M. A. and Kovacs, K. M. (2008),
#' A simple new algorithm to filter marine mammal Argos locations.
#' Marine Mammal Science, 24: 315?V325. doi: 10.1111/j.1748-7692.2007.00180.x

telem_sda <- function(data, ..., filter=TRUE){

  if( !is(data, "telemetry") & is(data, "list")){
    data <- lapply(data, telem_sda0, filter=filter,...)
  } else {
    data <- telem_sda0(data, ..., filter=filter)
  }

  return(data)
}

telem_sda0 <- function(data, ..., filter=TRUE){
  if(!check_telem(data)) stop("'data' does not appear to be a 'telemetry' object.")
  tr <- as(data[,c('longitude','latitude','timestamp')], "data.frame")
  tr$id <- 1
  tr <- suppressWarnings(trip(tr,  c("timestamp","id")))
  keep <- sda(tr,...)
  if(filter){
    data <- data[keep,]
  } else {
    data$keep <- keep
  }
  return(data)
}

