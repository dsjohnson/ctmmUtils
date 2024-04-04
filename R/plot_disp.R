#' @title Plot migration detection results
#' @description PLot the animal dispersion from the base location over time. Points
#' are colored to reflect estimated migration and non-migration phases.
#' See \code{\link{migration_det}}
#' @param data Original data used by `\link[ctmmUtils]{migration_det}` call.
#' @param migr_tbl Results table produced by `\link[ctmmUtils]{migration_det}`.
#' @param interactive Logical. If `TRUE` the `plotly` package will be used to make the plot more interactive.
#' @export
#' @author Devin S. Johnson
#' @import sf fuzzyjoin dplyr ggplot2
#' @importFrom plotly ggplotly
plot_disp <- function(data, migr_tbl, interactive=FALSE){
  individual.local.identifier <- NULL

  mult_anim <- TRUE
  if(!"individual.local.identifier" %in% colnames(data)){
    mult_anim <- FALSE
  } else{
    id <- unique(data$individual.local.identifier)
    if(length(id)==1) mult_anim <- FALSE
  }

  if(mult_anim){
    n <- length(id)
    plt <- vector("list", n)
    if(nrow(migr_tbl)!=n) stop("'migr_tbl' is not the same length as the number of unique individuals in 'data'")
    for(i in 1:n){
      data_id <-  filter(data, individual.local.identifier==id[[i]])
      plt[[i]] <- plot_disp0(data_id, migr_tbl$migr_tbl[[i]])
      plt[[i]] <- plt[[i]] + ggtitle(paste0("Animal: ",i))
      if(interactive) plt[[i]] <- ggplotly(plt[[i]])
    }
  } else{
    plt <- plot_disp0(data, migr_tbl)
    if(interactive) plt <- ggplotly(plt)
  }
  return(plt)
}


plot_disp0 <- function(data, migr_tbl){
  timestamp <- dist <- migr_evt <- NULL
  base <- attr(migr_tbl, "base")
  ddd <- data.frame(
    timestamp = data$timestamp,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  ddd <- fuzzyjoin::fuzzy_left_join(ddd,migr_tbl,
                                    by=c(timestamp="start",timestamp="end"),
                                    match_fun = list(`>=`, `<=`))
  plt <- ggplot() +
    geom_point(aes(x=timestamp, y=dist), alpha=1, color="slategray2",
               data=ddd %>% filter(migr_evt=="0")) +
    geom_point(aes(x=timestamp, y=dist), alpha=1, color="darkred",
               data=ddd %>% filter(migr_evt=="1")) +
    theme_classic() +
    xlab("Date") + ylab("Dispersal (km)")
  return(plt)
}
