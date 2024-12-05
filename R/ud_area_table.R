#' @title Get UD Area Table
#' @param uds A named list of `ctmm` `UD` objects
#' @param level.UD A vector of UD levels for computation
#' @param unit Units for area estimates. Defaults to `unit = "m^2"`
#' @export

ud_area_table <- function(uds, level.UD=0.95, unit=NULL){
  out <- NULL
  for(i in 1:length(uds)){
    for(j in 1:length(level.UD)){
      if(!inherits(uds[[i]], "UD")){
        tmp <- data.frame(individual.local.identifier = names(uds)[i], level_ud=level.UD[j], area=NA, low_ci=NA, hi_ci=NA)
        out <- rbind(out, tmp)
      } else{
        summ_u <- as.vector(summary(uds[[i]], convex=FALSE, level=0.95, level.UD=level.UD[j], units=FALSE)$CI)
        tmp <- data.frame(individual.local.identifier = names(uds)[i], level_ud=level.UD[j], area=summ_u[2], low_ci=summ_u[1], hi_ci=summ_u[3])
        out <- rbind(out, tmp)
      }
    }
  }
  units(out$area) <- as_units("m^2")
  units(out$low_ci) <- as_units("m^2")
  units(out$hi_ci) <- as_units("m^2")
  if(!is.null(unit)){
    units(out$area) <- as_units(unit)
    units(out$low_ci) <- as_units(unit)
    units(out$hi_ci) <- as_units(unit)
  }
  return(out)
}
