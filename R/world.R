#require(rnaturalearth)
#' @export

world <- function(){
  countries <- terra::vect('Q:\\Shared drives\\Data\\Original\\ne_10m_admin_0_countries_lakes\\ne_10m_admin_0_countries_lakes.shp')
  return(countries)
}
