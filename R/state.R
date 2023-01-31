#require(rnaturalearth)
#' @export

state <- function(name){
  usa <- rnaturalearth::ne_states(country='United States of America')
  ste <- usa[usa$name==name,]
  return(ste)
}
