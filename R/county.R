#require(rnaturalearth)
#' @export

county <- function(state, names){
  usa <- terra::vect('Q:\\Shared drives\\Data\\Vector\\USA\\us_lower_48_counties.gpkg')
  ste <- usa[usa$STATE_NAME==state,]
  cty <- ste[ste$NAME%in%names,]
  return(cty)
}
