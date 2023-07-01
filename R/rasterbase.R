#' @export

rasterbase <- function(res){
  geodir <- 'Q:\\Shared drives\\Data\\Raster\\'
  if(res==1000){
    if(file.exists(paste(geodir, 'Global\\base_1000m.tif', sep=''))==F){
      base.rast <- geodata::worldclim_global(var='bio', res=.5, path=paste(geodir, 'Global\\',sep=''))
      base.rast <- base.rast[[1]]; base.rast <- base.rast/base.rast
      terra::writeRaster(base.rast, paste(geodir, 'Global\\base_1000m.tif', sep=''))
    }
    if(file.exists(paste(geodir, 'Global\\base_1000m.tif', sep=''))==T){
      base.rast <- terra::rast(paste(geodir, 'Global\\base_1000m.tif', sep=''))
    }
  }
  if(res<1000){
    if(file.exists(paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))==F){
      if(res==500){
        base.rast <- terra::rast('Q:\\Shared drives\\Data\\Original\\BioClimComposite_1971_2000_400m.tif')
        base.rast <- base.rast[[1]]
        base.rast <- terra::project(base.rast, "epsg:4326", threads=T)
        base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
        base.rast <- terra::crop(base.rast, y=pops.sdm::l48(), mask=T)
      }
      if(res==250){
        base.rast <- terra::rast()
      }
      if(res==100){
        base.rast <- terra::rast(paste(geodir, 'Global\\pop\\GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif', sep=''))
        base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
        base.rast <- terra::crop(x=base.rast, y=pops.sdm::l48(), mask=T)
      }
      if(res==33){base.rast <- terra::rast()}
      base.rast <- base.rast/base.rast
      terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
    }
    if(file.exists(paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))==T){
      base.rast <- terra::rast(paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
    }
  }
  return(base.rast)
}
