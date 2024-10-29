#' @title Convert Wildlife Computers data imported with
#' `\link[ctmmUtils]{read_wc_dirs}` to a `telemetry` object from the
#' `ctmm` package.
#' @param x An sf data frame output by the function `\link[ctmmUtils]{read_wc_dirs}`.
#' @param ... Additional arguments to be passed to `\link[ctmm]{as.telemetry}`
#' @author Josh M. London, Devin S. Johnson
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom ctmm as.telemetry uere<-
#' @importFrom tidyr nest
#'
#' @export
as_telem <- function(x,...){
  type <- individual.local.identifier <- timestamp  <- NULL

  # browser()

  xargs <- list(...)
  if(!inherits(x, "sf")) stop("'x' not an {sf} data frame.")
  if(!st_is_longlat(x)){
    prj <- st_crs(x)$proj4string
    if(!is.null(xargs$projection)){
      message("\nAdditionally provided projection will override the current projection of the data.\n\n")
    } else{
      xargs$projection <- prj
    }
  } else{
    if(is.null(xargs$projection)) stop("A projection must be provided.")
  }

  x <- st_drop_geometry(x)

  mult_anim <- TRUE
  if(!"individual.local.identifier" %in% colnames(x)){
    mult_anim <- FALSE
  } else{
    id <- unique(x$individual.local.identifier)
    if(length(id)==1) mult_anim <- FALSE
  }

  if(!mult_anim){
    out <- as_telem0(x, xargs)
  } else{
    out <- tibble(individual.local.identifier = id, telem=vector("list",length(id)))
    for(i in 1:nrow(out)){
      locs_id <- filter(x, individual.local.identifier==id[[i]])
      out$telem[[i]] <- as_telem0(locs_id, xargs)
    }
  }
  return(out)
}

#' @importFrom ctmm as.telemetry tbind
as_telem0 <- function(x, xargs) {
  type <- NULL
  x$individual.local.identifier <- 0
  lv_quality <- levels(x$quality)
  x$quality <- as.character(x$quality)
  x <- rm_dup0(x)

  locs_f <- x |> filter(type%in%c("FastGPS","known"))
  locs_a <- x |> filter(type=="Argos")
  rm(x)

  if(nrow(locs_a)>0){
    locs_a <- locs_a |> st_drop_geometry()
    has_argos_mix <-  any(!is.na(locs_a$Argos.semi.major)) & any(is.na(locs_a$Argos.semi.major))
    if(has_argos_mix){
      locs_a_kf <- filter(locs_a, !is.na(locs_a$Argos.semi.major))
      locs_a_kf <- do.call(ctmm::as.telemetry, c(list(object=locs_a_kf), xargs))
      locs_a_dop <- filter(locs_a, is.na(locs_a$Argos.semi.major))
      locs_a_dop<- do.call(ctmm::as.telemetry, c(list(object=locs_a_dop), xargs))
      locs_a <- ctmm::tbind(locs_a_dop, locs_a_kf)
    } else{
      locs_a <- do.call(ctmm::as.telemetry, c(list(object=locs_a), xargs))
    }
  }
  if(nrow(locs_f)>0){
    locs_f <- locs_f %>% st_drop_geometry() %>% mutate(
      HDOP = dplyr::case_when(
        type == "known" ~ sqrt(2),
        type=="FastGPS" & quality=="4" ~ sqrt(2)*(1163)/20,
        type=="FastGPS" & quality=="5" ~ sqrt(2)*(169)/20,
        type=="FastGPS" & quality=="6" ~ sqrt(2)*(71)/20,
        type=="FastGPS" & quality=="7" ~ sqrt(2)*(43)/20,
        type=="FastGPS" & quality=="8" ~ sqrt(2)*(34)/20,
        type=="FastGPS" & quality=="9" ~ sqrt(2)*(28)/20,
        type=="FastGPS" & quality=="10" ~ sqrt(2)*(24)/20,
        type=="FastGPS" & quality=="11" ~ sqrt(2),
        TRUE ~ Inf
      )
    )
    locs_f <- do.call(ctmm::as.telemetry, c(list(object=locs_f), xargs))
    uere(locs_f) <- 20
  }

  if(nrow(locs_a)>0 & nrow(locs_f)>0) out <-suppressWarnings(tbind(locs_a, locs_f))
  if(nrow(locs_a)==0) out <- locs_f
  if(nrow(locs_f)==0) out <- locs_a

  if("quality"%in%colnames( out)){
    out$quality <- factor(out$quality, levels=lv_quality)
  }

  return(out)
}
