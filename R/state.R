#require(rnaturalearth)
#' @export

state <- function(names){
  usa <- rnaturalearth::ne_states(country='United States of America')
  ste <- usa[usa$name%in%names,]
  return(terra::vect(ste))
}
