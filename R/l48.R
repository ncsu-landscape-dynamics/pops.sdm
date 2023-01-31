#require(rnaturalearth)
#' @export

l48 <- function(){
  bbox <- c(24.5, -125, 49.5, -66.5)
  us_can <- borders[borders$SOVEREIGNT%in%c('United States of America', 'Canada'),]
  us_can <- terra::crop(us_can, ext(c(bbox[2], bbox[4], bbox[1], bbox[3])))
  return(us_can)
}
