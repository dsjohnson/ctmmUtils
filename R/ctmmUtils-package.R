#' @title Functions To Increase Usability Of The \code{ctmm} Package
#'
#' @description This package is a collection of functions that enhance the \code{ctmm} package for
#' for analysis of animal telemetry data.
#'
#' \tabular{ll}{
#' Package: \tab ctmmUtils\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.0.9106\cr
#' Date: \tab December 3, 2024\cr
#' License: \tab CC0 \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @note This software package is developed and maintained by scientists at the
#' NOAA Fisheries Pacific Islands Fisheries Science Center and should be
#' considered a fundamental research communication. The recommendations and
#' conclusions presented here are those of the authors and this software should
#' not be construed as official communication by NMFS, NOAA, or the U.S. Dept.
#' of Commerce. In addition, reference to trade names does not imply endorsement
#' by the National Marine Fisheries Service, NOAA. While the best efforts have
#' been made to insure the highest quality, tools such as this are under
#' constant development and are subject to change.
#'
#' @name ctmmUtils-package
#' @aliases ctmmUtils
#' @author Devin S. Johnson and Josh M. London
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#'
"_PACKAGE"



.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(package, version, paste("(",date, ")", sep=""))
  )
}


#' @importFrom methods is
check_telem <- function(data){
  chk <- FALSE
  if(!is(data, "telemetry") & is(data,"list")){
    chk <- all(sapply(data, is,  "telemetry"))
  } else{
    chk <- is(data,"telemetry")
  }
  return(chk)
}

#' @import dplyr
rm_dup <- function(x){
  individual.local.identifier <- timestamp <- quality <- NULL
  x <- x |>
    group_by(individual.local.identifier) |>
    arrange(timestamp, quality) |>
    mutate(
      rank = 1L,
      rank = case_when(
        duplicated(timestamp, fromLast = FALSE) ~ lag(rank) + 1L,
        TRUE ~ rank)
      ) |>
    dplyr::filter(rank == 1) |>
    ungroup() |>
    arrange(individual.local.identifier, timestamp)
}

#' @import dplyr
rm_dup0 <- function(x){
  timestamp <- quality <- NULL
  x <- x |>
    arrange(timestamp, quality) |>
    mutate(
      rank = 1L,
      rank = case_when(
        duplicated(timestamp, fromLast = FALSE) ~ lag(rank) + 1L,
        TRUE ~ rank)
    ) |>
    dplyr::filter(rank == 1) |>
    arrange(timestamp)
}


# rename_geometry <- function(g, name){
#   current = attr(g, "sf_column")
#   names(g)[names(g)==current] = name
#   st_geometry(g)=name
#   g
# }


# weighted.var <- function(x, w, na.rm = FALSE) {
#   if (na.rm) {
#     w <- w[i <- !is.na(x)]
#     x <- x[i]
#   }
#   sum.w <- sum(w)
#   sum.w2 <- sum(w^2)
#   mean.w <- sum(x * w) / sum(w)
#   (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
#                                        na.rm)
# }

#~~~ Reference for summarizing overlapping polys:
# https://stackoverflow.com/questions/48279545/summarise-attributes-from-sfst-intersection-where-geometries-overlaps
# *** This doesn't really work in practice. Just use a common grid. ***






