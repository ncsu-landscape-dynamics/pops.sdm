geodir <- 'Q:\\Shared drives\\Data\\Raster\\'


### climate 400m
if(bio=T){
bioc.400 <- terra::rast(paste(geodir, 'Global\\BioClimComposite_1971_2000_400m.tif', sep=''))
}

#### soil 250m
if(soil=T){
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
if(pop=T){
popd <- terra::rast()
}

#  elev1m <-
