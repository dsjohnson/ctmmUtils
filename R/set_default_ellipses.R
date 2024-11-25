#' @title Set Default Values for MIssing ARGOS Diagnostic Data
#' @param x An sf data frame output by the function `\link[ctmmUtils]{read_wc_dirs}`.
#' @description
#' This function will fill in missing ARGOS diagnostic ellipse data for either type of location, `Argos` or `FastGPS`. This is necessary for converting the
#' data to `telemetry` format for the `ctmm` package. The default values used for non-KF Argos data (old location classes) are those in Vincent et al. (2002; Table 1 nonfilterd). For
#' `FastGPS` data, the values in Dujon et al. (2014; Table 1, 95\%-tile unfiltered.).
#'
#' Because Dujon et al. (2014) only used error magnitude, we assumed normally distributed error positions, \eqn{N(0,\sigma^2)}. Thus, the error magitudes are Rayleigh distributed. So
#' we used the CDF of the Rayleigh distribution and the 95\%-tile given in the table to solve for the appropriate \eqn{\sigma}.
#'
#'@importFrom dplyr mutate case_when bind_rows
#'@export

set_default_ellipses <- function(x){

type <- Argos.semi.major <- Argos.semi.minor <- Argos.location.class <- timestamp <- NULL

N <- nrow(x)
miss_idx <- with(x, is.na(Argos.semi.major) | is.na( Argos.semi.minor) | is.na(Argos.orientation))
miss_diag <- x[miss_idx,]
has_diag <- x[!miss_idx,]

miss_diag <- miss_diag |> mutate(
  Argos.semi.major = dplyr::case_when(
    type == "known" ~ sqrt(2),
    !is.na(Argos.error.radius) ~ Argos.error.radius,
    quality=="4" ~ sqrt(2)*475.1308,
    quality=="5" ~ sqrt(2)*69.04309,
    quality=="6" ~ sqrt(2)*29.00627,
    quality=="7" ~ sqrt(2)*17.56718,
    quality=="8" ~ sqrt(2)*13.89033,
    quality=="9" ~ sqrt(2)*11.43909,
    quality=="10" ~ sqrt(2)*9.804936,
    quality=="11" ~ sqrt(2)*7.762241,
    quality=="B"  ~ sqrt(2)*41219/2,
    quality=="A"  ~ sqrt(2)*10393/2,
    quality=="0"  ~ sqrt(2)*15361/2,
    quality=="1"  ~ sqrt(2)*3498/2,
    quality=="2"  ~ sqrt(2)*1355/2,
    quality=="3"  ~ sqrt(2)*742/2,
    TRUE ~ NA
  ),
  Argos.semi.minor = dplyr::case_when(
    type == "known" ~ sqrt(2),
    !is.na(Argos.error.radius) ~ Argos.error.radius,
    quality=="4" ~ sqrt(2)*475.1308,
    quality=="5" ~ sqrt(2)*69.04309,
    quality=="6" ~ sqrt(2)*29.00627,
    quality=="7" ~ sqrt(2)*17.56718,
    quality=="8" ~ sqrt(2)*13.89033,
    quality=="9" ~ sqrt(2)*11.43909,
    quality=="10" ~ sqrt(2)*9.804936,
    quality=="11" ~ sqrt(2)*7.762241,
    quality=="B" ~ sqrt(2)*15535/2,
    quality=="A" ~ sqrt(2)*5373/2,
    quality=="0" ~ sqrt(2)*5517/2,
    quality=="1" ~ sqrt(2)*1265/2,
    quality=="2" ~ sqrt(2)*511/2,
    quality=="3" ~ sqrt(2)*326/2,
    TRUE ~ NA
  ),
  Argos.orientation = ifelse(type=="Argos", 90, 0),
  Argos.error.radius = Argos.semi.major,
  Argos.location.class = ifelse(type=="FastGPS", 3, Argos.location.class)
)

out <- bind_rows(miss_diag, has_diag) %>% arrange(timestamp)
if(nrow(out)!=N) warning("Some locations seem to have been lost please compare output to the original data.")
if(any(is.na(out$Argos.semi.major)) | any(is.na(out$Argos.semi.minor)) | any(is.na(out$Argos.orientation))) warning("There are still observations with missing Argos ellipse data, please check!")
return(out)
}
