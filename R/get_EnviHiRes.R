# require(geodata)
# require(terra)
# require(parallel)
#' @export

get_EnviHiRes <- function(bio=F, elev=F, gdd=F, lc=F, pop=F, ptime=F, rnr=F, soil=F, tbase=5){


geodir <- 'Q:\\Shared drives\\Data\\Raster\\'


### climate 400m
if(bio==T){
biovar <- terra::rast(paste(geodir, 'Global\\BioClimComposite_1971_2000_400m.tif', sep=''))
}

#### soil 250m
if(soil==T){
if(grepl('', list.files(paste(geodir, 'Global', sep='')))){
soil.h20.1500 <- terra::rast(paste(geodir, 'Global\\soils\\250m\\1500kpa.mean.tif', sep=''))
soil.h20.33 <- terra::rast(paste(geodir, 'Global\\soils\\250m\\33kpa.mean.tif', sep=''))
soil.ph <-  terra::rast(paste(geodir, 'Global\\soils\\250m\\ph.mean.tif', sep=''))
}
# type <- c('ph.h2o', '1500kPa', '33kPa')
# i <- 1
# for(i in i:length(type)){
#   sl.list <- lapply(list.files(path='C:\\Users\\bjselige\\Desktop\\Soil_OpenLandMap\\', pattern=type[i], full.names=T), raster)
#   sl.mean <- terra::mean(sl.list[[1]], sl.list[[2]], sl.list[[3]], sl.list[[4]], sl.list[[5]], sl.list[[6]])
#   writeRaster(sl.mean, 'C:\\Users\\bjselige\\Desktop\\Soil_OpenLandMap\\33kPa.mean.tif')
# }
}

#### land cover 30m
# if(land=T){
# nlcd.30 <- terra::rast(paste(geodir, , sep=''))
# impsurf <- terra::rast(paste(geodir, , sep=''))
# }

### population density
if(pop==T){
popd <- terra::rast()
}

#  elev1m <-

if(bio==F){biovar <- NULL; biocl <- NULL}
if(elev==F){elevvar <- NULL; elevcl <- NULL}
if(gdd==F){gddvar <- NULL; gddcl <- NULL}
if(lc==F){lcvar <- NULL; lccl <- NULL}
if(pop==F){popvar <- NULL; popcl <- NULL}
if(ptime==F){timevar <- NULL; timecl <- NULL}
if(rnr==F){rnrvar <- NULL; rnrcl <- NULL}
if(soil==F){solvar <- NULL; solcl <- NULL}

envi <- terra::rast(terra::as.list(c(biovar, elevvar, gddvar, lcvar, popvar, timevar, rnrvar, solvar)))
clst <- rbind(biocl, elevcl, gddcl, lccl, popcl, timecl, rnrcl, solcl)
clst$cluster <- as.integer(as.factor(clst$cluster))
cl2 <- clst$cluster; names(cl2) <- gsub(' ', '.', clst$var)
return(list('rast'=envi, 'clust'=cl2))

}
