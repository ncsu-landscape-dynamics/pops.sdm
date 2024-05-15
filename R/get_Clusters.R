# require(klaR)
# require(terra)
#' @export

get_Clusters <- function(envi, n=NULL, samp=NULL){
  if(!is.null(samp)){
    # ev.nonna <- !is.na(envi[[1]])
    # ev.cells <- which(terra::values(ev.nonna)==T)
    # ev.samp <- sample(x=ev.cells, size=100000000)
    if(samp>90000000){samp <- 90000000}
    envi.x <- exactextractr::exact_extract(envi, as(domain, 'Spatial'), force_df=T)
    envi.b <- dplyr::bind_rows(envi.x); envi.b <- envi.b[,1:(length(envi.b)-1)]
    remove(envi.x)
    envi.df <- na.omit(envi.b)
    envi.df2 <- envi.df[sample(1:nrow(envi.df), size=samp),]
  }
  if(is.null(samp)){
    envi.df2 <- terra::values(envi, na.rm=T, dataframe=T)
  }
  envi.cc <- klaR::corclust(x=envi.df2, method='complete') #envi.pl <- plot(envi.cc, mincor=.7, selection='numeric')
  if(!is.null(n)){envi.tr <- klaR::cvtree(envi.cc, k=n)}
  if(is.null(n)){envi.tr <- klaR::cvtree(envi.cc, mincor=.1)}
  return(envi.tr)
}
