#' @title Convert Wildlife Computers data imported with
#' `\link[ctmmUtils]{read_wc_dirs}` to a `telemetry` object from the
#' `ctmm` package.
#' @param x An sf data frame output by the function `\link[ctmmUtils]{read_wc_dirs}`.
#' @param ... Additional arguments to be passed to `\link[ctmm]{as.telemetry}`
#' @author Josh M. London, Devin S. Johnson
#' @import dplyr
#' @importFrom ctmm as.telemetry
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

  if("individual.local.identifier" %in% colnames(x)){
    x <- split(x, x$individual.local.identifier)
    tdata <- lapply(x, as_telem0, xargs=xargs)
    names(tdata) <- names(x)
  } else{
    tdata <- as_telem0(x,xargs)
  }

  if(inherits(tdata, "list") && length(tdata)==1) tdata <- tdata[[1]]

  return(tdata)
}





as_telem0 <- function(x0, xargs){
  covs <- argosDiag2Cov(x0$Argos.semi.major, x0$Argos.semi.minor, x0$Argos.orientation)
  covs <- zapsmall(covs)
  tdata0 <- do.call(ctmm::as.telemetry, c(list(object=x0), xargs))
  tdata0$COV.x.x <- covs$Cov.x.x
  tdata0$COV.x.y <- covs$Cov.x.y
  tdata0$COV.y.y <- covs$Cov.y.y
  return(tdata0)
}


argosDiag2Cov <- function(Major, Minor, Orientation){
  a <- Major
  b <- Minor
  if(any(b<=.Machine$double.eps, na.rm=TRUE)) stop("There are very small (or 0) values for the minor ellipse lengths! These may need to be removed.")
  theta <- Orientation
  if(any(theta < 0 | theta > 180, na.rm = TRUE)) stop("Argos diagnostic data orientation outside of [0,180]!")
  if(any(a < 0, na.rm = TRUE)) stop("Argos diagnostic data major axis < 0!")
  if(any(b < 0, na.rm = TRUE)) stop("Argos diagnostic data minor axis < 0!")
  theta <- pi*(theta/180)
  k <- sqrt(2)
  v1 <- (a/k)^2*sin(theta)^2 + (b/k)^2*cos(theta)^2
  v2 <- (a/k)^2*cos(theta)^2 + (b/k)^2*sin(theta)^2
  c12 <- ((a^2 - b^2)/k^2)*cos(theta)*sin(theta)
  return(data.frame(Cov.x.x=v1, Cov.x.y=c12, Cov.y.y=v2))
}





