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
    envi.b <- dplyr::bind_rows(envi.x); envi.b <- envi.b[,1:(length(envi.b)-1)] #removes 'coverage fraction'
    remove(envi.x)
    envi.d <- na.omit(envi.b)
    remove(envi.b)
    envi.df <- envi.d[sample(1:nrow(envi.d), size=samp),]
    if(any(grepl('landcover', names(envi.df)))){
      envi.df$landcover <- factor(envi.df$landcover,
                                  levels=c(12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
                                  labels=c('Ice', 'Dev_1', 'Dev_2', 'Dev_3', 'Dev_4', 'Barren', 'Decid', 'Everg',
                                    'Mixed', 'Shrub', 'Grass', 'Pastr', 'Culti', 'WetWdy', 'WetHrb'))
    }
  }
  if(is.null(samp)){
    envi.df <- terra::values(envi, na.rm=T, dataframe=T)
  }
  envi.cc <- klaR::corclust(x=envi.df) #envi.pl <- plot(envi.cc, mincor=.7, selection='numeric')
  if(!is.null(n)){envi.tr <- klaR::cvtree(envi.cc, k=n)}
  if(is.null(n)){envi.tr <- klaR::cvtree(envi.cc, mincor=.1)}
  return(envi.tr)
}
