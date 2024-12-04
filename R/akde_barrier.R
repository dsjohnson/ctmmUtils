#' @title Safe method for calculating an AKDE with a barrier
#' @param data 2D timeseries telemetry data represented as a telemetry object or list of objects.
#' @param CTMM A ctmm movement model from the output of ctmm.fit or list of objects.
#' @param R A named list of raster covariates if CTMM contains an RSF model.
#' @param VMM An optional vertical ctmm object for 3D home-range calculation.
#' @param barrier An `sf` polygon data frame for enforcing hard boundaries
#' @param barrier_in Locations are assumed to be inside the `sf` polygons if `barrier_in=TRUE`
#'  and outside of `barrier` if `barrier_in=FALSE`.
#' @param debias Debias the distribution for area estimation (AKDEc).
#' @param smooth "Smooth" out errors from the data.
#' @param weights Optimally weight the data to account for sampling bias (See `\link[ctmm]{bandwidth}` for akde details).
#' @param error Target probability error.
#' @param res Number of grid points along each axis, relative to the bandwidth.
#' @param grid Optional grid specification via raster, UD, or list of arguments (See ‘Details’ in `\link[ctmm]{akde}`).
#' @param ... Arguments passed to akde, `\link[ctmm]{bandwidth}`, and `\link[ctmm]{mean.ctmm}`.
#' @importFrom ctmm raster akde
#' @importFrom terra rast mask as.matrix vect values `values<-`
#' @importFrom sf st_transform st_crs st_union
#' @export

akde_barrier <- function(data, CTMM, VMM=NULL, R=list(), barrier=NULL, barrier_in=FALSE,
                              debias=TRUE, weights=FALSE, smooth=TRUE, error=0.001,
                              res=10,grid=NULL,...){
  ud <- akde(data=data, CTMM=CTMM, VMM=VMM, R=R, debias=debias, weights=weights, smooth=smooth, error=error, res=res, grid=grid,...)
  nr <- length(ud$r$y)
  a <-  prod(ud$dr)
  r_cdf <- terra::rast(raster(ud, "CDF"))
  r_pdf <- terra::rast(raster(ud, "PDF"))
  barrier <- st_transform(barrier, st_crs(r_cdf)) |> st_union()
  r_cdf <- mask(r_cdf, vect(barrier), inverse=!barrier_in, updatevalue=1)
  r_pdf <- mask(r_pdf, vect(barrier), inverse=!barrier_in, updatevalue=0)
  values(r_pdf) <- values(r_pdf)/(sum(values(r_pdf)) * a)
  ud$CDF <-  t(as.matrix(r_cdf, wide=TRUE)[nr:1,])
  ud$PDF <-  t(as.matrix(r_pdf, wide=TRUE)[nr:1,])
  return(ud)
}
