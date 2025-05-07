#' @title Read individual telemetry data from Wildlife Computers portal directories
#' @description Read and combine data downloaded from Wildlife Computers portal
#' into individual directories.
#' @param x Directory containing the individual telemetry data directories.
#' @param remove_duplicates Logical. Should observations with duplicated times be removed? The observation
#' with the highest quality will be retained.
#' @param default_ellipse Logical. Should a default set of Argos error ellipse values be added for `FastGPS` and `Argos` doppler locations?
#' See `\link[ctmmUtils]{set_default_ellipses}` for further information.
#' @param rm_mote Logical, default = `TRUE`. If set to `rm_mote = TRUE` the locations that are `type == "Mote"` stations will be
#' removed, due to lack of location error information. To force retension, set `rm_mote = FALSE`.
#' @export
#' @author Devin S. Johnson, Josh M. London
#' @import lubridate dplyr sf
#' @importFrom janitor clean_names
#' @importFrom readr read_csv

read_wc_dirs <- function(x, remove_duplicates=TRUE, default_ellipse=TRUE, rm_mote=TRUE){

  Longitude <- Latitude <- quality <- deploy_id <- datetime <- longitude <- latitude <-
    error_ellipse_orientation <- error_semi_minor_axis <- error_semi_major_axis <-
    type <- location.long <- location.lat <- error_radius <- NULL
  # Determine which file to load for each animal:
  dirs <- list.dirs(x)[-1]

  # Read in data and combine into single table
  # There are 2 animals with no location data
  locs <- NULL
  no_dat <- NULL
  for(i in 1:length(dirs)){
    loc_paths <- list.files(dirs[[i]], pattern="*-Locations.csv")
    if(length(loc_paths)==0){
      no_dat <- c(no_dat, dirs[[i]])
      next
    }
    if(length(loc_paths)==1){
      loc_file <- loc_paths[[1]]
    } else{
      loc_paths <- strsplit(loc_paths,"-")
      loc_paths <- loc_paths[sapply(loc_paths, "length")>2]
      run <- as.numeric(sapply(loc_paths, \(x) x[[2]]))
      loc_file <- paste(loc_paths[run==max(run)][[1]], collapse="-")
    }
    loc_file <- paste0(dirs[[i]],"/",loc_file)
    id_data <- readr::read_csv(loc_file, show_col_types=FALSE) %>%
      filter(!is.na(Latitude), !is.na(Longitude))
    time <- parse_date_time(id_data$Date,"%H:%M:%S %d-%b-%Y", quiet=TRUE)
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%y %H:%M", quiet=TRUE)
    }
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%y %H:%M:%S", quiet=TRUE)
    }
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%Y %H:%M:%S", quiet=TRUE)
    }
    if(all(is.na(time))){
      stop("Date format is unrecognizable.")
    }
    if(any(is.na(time))){
      warning("Some dates were not converted to POSIX format. Look for NAs in datetime.")
    }
    id_data$datetime <- time
    id_data <- select(id_data, -Date)
    id_data <- janitor::clean_names(id_data)



    locs <- bind_rows(locs, id_data)
  }

  locs <- locs |> mutate(
    quality = ifelse(type=="Mote", "M", quality),
    quality = factor(quality, levels=c(as.character(11:0),"A","B","Z","M")),
    high_qual = ifelse(!quality %in% c("4","0","A","B","Z","M"), 1, 0)
  )

  # rename for as.telemetry
  locs <- locs |>
    rename(
      individual.local.identifier = deploy_id,
      timestamp = datetime,
      location.long = longitude,
      location.lat = latitude,
      Argos.error.radius = error_radius,
      Argos.orientation = error_ellipse_orientation,
      Argos.semi.minor = error_semi_minor_axis,
      Argos.semi.major = error_semi_major_axis
    ) %>% mutate(
      Argos.location.class = ifelse(type=="Argos", as.character(quality), NA),
      x = location.long,
      y = location.lat
    )

  if(any(locs$type=="Mote") & rm_mote) locs <- locs[locs$type!="Mote",]
  if(remove_duplicates) locs <- rm_dup(locs)
  locs <- locs[locs$quality!="Z",]
  if(default_ellipse) locs <- set_default_ellipses(locs)

  locs <- droplevels(locs)

  # if(any(is.na(locs$x)) | any(is.na(locs$y))) browser()
  locs <- st_as_sf(locs, coords = c("x","y"), crs=4326)

  if(!is.null(no_dat)){
    www <- paste0("There are individuals without location data files: \n", paste(no_dat, "\n", collapse = ""))
    warning(www, immediate. = TRUE)
  }

  return(locs)
}



