#require(rnaturalearth)
#' @export

country <- function(names){
  country <- rnaturalearth::ne_countries(scale=10, country=names)
  return(terra::vect(country))
}
