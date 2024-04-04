#' @title Migration detection
#' @description Creates a data table that indicates the times of different
#' phases of movement. This method uses changes in the overall dispersion rate of the animal
#' from the 'base' time to detect changes in overall movement from small scale local movement
#' to large scale migration.
#' @param data A `telemytry` data object of locations (see `\link[ctmm]{as.telemetry}`)
#' @param min_disp The minimum dispersion rate to be considered a migration interval,
#' e.g. 10 for a 10km dispersion minimum.
#' @param max_num_mig The maximum number of migration intervals.
#' @param min_phase_len The minimum length of time that a migration or non-migration event
#' will take, e.g., 7 implies a minimum of 7 time intervals for a phase.
#' @param grid_res The temporal resolution at which migrations are detected. e.g., "day" (default) implies migration
#' start and end is detected on a daily resolution.
#' @param base The location at which dispersion is measured. Can be one of \code{"first"} (first location),
#' \code{"last"} (final location), or some other \code{sf::sfc} point location.
#' @param max_k The maximum degrees of freedom used by \code{mgcv::gam} to model dispersion and estimate
#' the derivative of the dispersion function.
#' @export
#' @author Devin S. Johnson
#' @import units lubridate dplyr sf
#' @rawNamespace import(mgcv, except = mvn)
#' @importFrom stats coef dist predict vcov
#' @importFrom methods is
#' @importFrom ctmm as.sf
migration_det <- function(data, min_disp, max_num_mig=1, min_phase_len=3,
                          grid_res="day", base="first", max_k=100){
  individual.local.identifier <- NULL

  mult_anim <- TRUE
  if(!"individual.local.identifier" %in% colnames(data)){
    mult_anim <- FALSE
  } else{
    id <- unique(data$individual.local.identifier)
    if(length(id)==1) mult_anim <- FALSE
  }

  if(!mult_anim){
    out <- migration_det0(data, min_disp=min_disp, max_num_mig=max_num_mig,
                          min_phase_len=min_phase_len, grid_res=grid_res,
                          base=base, max_k=max_k)
    return(out)
  } else {
    out <- tibble(individual.local.identifier = id, migr_tbl=vector("list",length(id)))
    for(i in 1:nrow(out)){
      data_id <- filter(data, individual.local.identifier==id[[i]])
      out$migr_tbl[[i]] <- migration_det0(data_id, min_disp=min_disp, max_num_mig=max_num_mig,
                                          min_phase_len=min_phase_len, grid_res=grid_res,
                                          base=base, max_k=max_k)
    }
    return(out)
  }
}



#' @importFrom utils tail
migration_det0 <- function(data, min_disp, max_num_mig=1, min_phase_len=3,
                           grid_res="day",
                           base="first", max_k=100){
  migr_evt <- phase <- timestamp <- disp_rate <- NULL
  if(base=="first"){
    base <- data[1,] %>% st_geometry()
  } else if(base=="last"){
    base <- data[nrow(data),] %>% st_geometry()
  } else if(!inherits(base,c("sf","sfc"))){
    stop("Unrecognized 'base' argument!")
  }
  if(nrow(data)<= (max_k+1) ){
    k <- floor(nrow(data)/2)
  } else{
    k <- as.integer(max_k)
  }
  ddd <- data.frame(
    timestamp = data$timestamp,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  ddd$time <- with(ddd,
                   as.numeric(timestamp)/as.numeric(duration(1, grid_res))
  )
  suppressWarnings(dfit <- mgcv::gam(dist~s(time, k=k, bs='ad'), data=ddd, method='REML'))
  grid1 <- floor_date(min(ddd$timestamp), grid_res)
  grid2 <- ceiling_date(max(ddd$timestamp), grid_res)
  newdata <- data.frame(timestamp = seq(grid1, grid2, grid_res))
  newdata$time <- with(newdata,
                       as.numeric(timestamp)/as.numeric(duration(1, grid_res))
  )
  X <- predict(dfit, newdata=newdata, type="lpmatrix")
  der <- diff(X%*%coef(dfit))
  D <- -1*diag(nrow(newdata))
  for(i in 1:(nrow(D)-1)){
    D[i,i+1] <- 1
  }
  D <- D[-nrow(D),]
  V <- sqrt(diag(D%*%X%*%vcov(dfit, unconditional=T)%*%t(X)%*%t(D)))
  derup <- der + 1.96*V
  derlo <- der - 1.96*V
  sig <- (derup*derlo >0)*(abs(der)>=min_disp)

  x <- rle(as.numeric(sig))
  x$values[x$lengths<=min_phase_len & x$values==1] <- 0
  x$values[x$lengths<=min_phase_len & x$values==0] <- 1
  x <- rle(inverse.rle(x))
  top_runs <- sort(x$lengths[x$values==1], decreasing=TRUE)[1:max_num_mig]
  trl <- length(x$lengths[x$values==1 & x$lengths%in%top_runs])
  if(trl>max_num_mig) warning("There were multiple migration events with the same length")
  x$values[x$values==1 & x$length<min(top_runs)] <- 0
  x$values <- cumsum(x$values)*x$values
  newdata$migr_evt <- c(inverse.rle(x),tail(x$values,1))
  newdata$disp_rate <- c(der, NA)
  newdata$phase <-with(rle(newdata$migr_evt), rep(seq_along(values), lengths))
  summ <- group_by(newdata, migr_evt, phase) %>%
    summarize(
      start = min(timestamp),
      end = max(timestamp),
      avg_disp_rate = mean(disp_rate, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(phase)
  summ$end[1:(nrow(summ)-1)] <- summ$end[1:(nrow(summ)-1)]+duration(grid_res)
  summ$avg_disp_rate <- units::set_units(summ$avg_disp_rate, paste0("km/",grid_res), mode='standard')
  summ$migr_evt <- factor(summ$migr_evt)
  summ$phase <- factor(summ$phase)

  attr(summ, "base") <- base

  # class(summ) <- c(class(summ),"migr_tbl")

  return(summ)
}
