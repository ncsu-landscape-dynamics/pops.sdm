#require(rnaturalearth)
#' @export

country <- function(name){
  country <- rnaturalearth::ne_countries(scale=10, country=name)
  return(country)
}
