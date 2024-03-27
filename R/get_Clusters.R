# require(klaR)
# require(terra)
#' @export

get_Clusters <- function(envi, n=NULL){
  #ev.samp <- sample.int(n=terra::ncell(envi[[1]]), size=100000)
  #envi.df <- terra::extract(envi, ev.samp)
  envi.df <- terra::values(envi, na.rm=T, dataframe=T)
  envi.cc <- klaR::corclust(x=envi.df)
  #envi.pl <- plot(envi.cc, mincor=.7, selection='numeric')
  if(!is.null(n)){envi.tr <- klaR::cvtree(envi.cc, k=n)}
  if(is.null(n)){envi.tr <- klaR::cvtree(envi.cc, mincor=.1)}
  return(envi.tr)
}
