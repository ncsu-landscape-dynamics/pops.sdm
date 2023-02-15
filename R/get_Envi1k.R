# require(geodata)
# require(terra)
# require(parallel)
#' @export

get_Envi1k <- function(bio=F, elev=F, gdd=F, lc=F, pop=F, ptime=F, rnr=F, soil=F, tbase=5){

  geodir <- 'Q:\\Shared drives\\Data\\Raster\\Global\\'

  if(bio==T){
    biovar <- geodata::worldclim_global(var='bio', res=.5, path=geodir)
    names(biovar) <- c('Mean.Annual.Temp', 'Mean.Diurnal.Range', 'Isothermality', 'Temp.Seasonality',
                       'Max.Temp.Warmest.Month', 'Min.Temp.Coldest.Month', 'Temp.Annual.Range',
                       'Mean.Temp.Wettest.Quarter', 'Mean.Temp.Driest.Quarter', 'Mean.Temp.Warmest.Quarter',
                       'Mean.Temp.Coldest.Quarter', 'Annual.Precip', 'Precip.Wettest.Month',
                       'Precip.Driest.Month', 'Precip.Seasonality', 'Precip.Wettest.Quarter',
                       'Precip.Driest.Quarter', 'Precip.Warmest.Quarter', 'Precip.Coldest.Quarter')
    biocl <- data.frame(var=names(biovar), cluster=NA)
    biocl$cluster[which(biocl$var%in%c('Mean.Annual.Temp',
                                       'Max.Temp.Warmest.Month',
                                       'Mean.Temp.Wettest.Quarter',
                                       'Mean.Temp.Warmest.Quarter'))] <- 'Temp 1'
    biocl$cluster[which(biocl$var%in%c('Mean.Diurnal.Range',
                                       'Isothermality',
                                       'Temp.Seasonality',
                                       'Min.Temp.Coldest.Month',
                                       'Temp.Annual.Range',
                                       'Mean.Temp.Driest.Quarter',
                                       'Mean.Temp.Coldest.Quarter'))] <- 'Temp 2'
    biocl$cluster[which(biocl$var%in%c('Annual.Precip',
                                       'Precip.Wettest.Month',
                                       'Precip.Wettest.Quarter',
                                       'Precip.Warmest.Quarter',
                                       'Precip.Coldest.Quarter'))] <- 'Precip 1'
    biocl$cluster[which(biocl$var%in%c('Precip.Seasonality',
                                       'Precip.Driest.Month',
                                       'Precip.Driest.Quarter'))] <- 'Precip 2'
  }

  if(elev==T){
    elevvar <- geodata::elevation_global(res=.5, path=geodir); names(elevvar) <- 'Elevation'
    elevcl <- data.frame(var=names(elevvar), cluster='Elevation 1')
  }

  if(gdd==T){
    if(exists('tbase')==F){tbase <- 5}
    getGDD <- function(tbase){
      gdpath <- paste(geodir, 'gdd.base', tbase, '.tif', sep='')
      if(file.exists(gdpath)){g1 <- terra::rast(gdpath)}
      if(!file.exists(gdpath)){print('Calculating GDD')
        tavg <- geodata::worldclim_global(var='tavg', res=.5, path=geodir)
        g1 <- terra::app(x=tavg, tbase, fun=function(x, tbase){
          days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
          xbase <- x-tbase; xbase[which(xbase<0)] <- 0
          return(sum(xbase*days))}, cores=detectCores()/2)
        terra::writeRaster(g1, filename=gdpath); closeAllConnections()
      }
      names(g1) <- paste('GDD ', tbase, sep=''); return(g1)
    }
    gddvar <- getGDD(tbase)
    gddcl <- data.frame(var=names(gddvar), cluster='Temp 1')
  }

  if(lc==T){
    globe <- geodata::worldclim_global(var='bio', res=.5, path=geodir)
    built <- terra::extend(x=landcover(var='built', path=geodir), y=ext(globe))
    cropl <- terra::extend(x=landcover(var="cropland", path=geodir), y=ext(globe))
    grass <- terra::extend(x=landcover(var='grassland', path=geodir), y=ext(globe))
    shrub <- terra::extend(x=landcover(var='shrubs', path=geodir), y=ext(globe))
    trees <- terra::extend(x=landcover(var='trees', path=geodir), y=ext(globe))
    wetld <- terra::extend(x=landcover(var='wetland', path=geodir), y=ext(globe))
    lcvar <- c(built, cropl, grass, shrub, trees, wetld)
    names(lcvar) <- c('Built', 'Cropland', 'Grassland', 'Shrubs', 'Trees', 'Wetland')
    lccl <- data.frame(var=names(lcvar), cluster='Landcover 1')
  }

  if(pop==T){
    popvar <- geodata::population(year='2020', res=.5, path=geodir)
    names(popvar) <- 'Population'
    popcl <- data.frame(var=names(popvar), cluster='Population 1')
  }

  if(ptime==T){
    getPrecipTiming <- function(){
      ptpath <- paste(geodir, '\\precip.timing.tif', sep='')
      if(file.exists(ptpath)){prect <- rast(ptpath)}
      if(!file.exists(ptpath)){print('Calculating Precip Timing (DJF-JJA)')
        prec <- geodata::worldclim_global(var='prec', res=.5, path=geodir)
        prect <- terra::app(x=prec, cores=detectCores()/2,
                            fun=function(x){return(sum(x[[12]], x[[1]], x[[2]])-sum(x[[6]], x[[7]], x[[8]]))})
        terra::writeRaster(prect, filename=ptpath); closeAllConnections()
      }
      names(prect) <- 'Precip.Timing'; return(prect)
    }
    timevar <- getPrecipTiming()
    timecl <- data.frame(var=names(timevar), cluster='Precip 2')
  }

  if(rnr==T){
    rnrvar <- terra::rast(list.files(geodir, '.dist.wrld', full.names = T))
    names(rnrvar) <- c('Road.Dist', 'Rail.Dist')
    rnrcl <- data.frame(var=names(rnrvar), cluster='Roads/Rails 1')
  }

  if(soil==T){
    soil.files <- c('Soil_pH_0cm.tif', 'Soil_pH_mean.tif', 'Soil_pH_200cm.tif',
                    'Soil_h2o_33kpa.0cm.tif', 'Soil_h2o_33kpa_mean.tif',
                    'Soil_h2o_33kpa_200cm.tif', 'Soil_h2o_1500kpa_0cm.tif',
                    'Soil_h2o_1500kpa_mean.tif', 'Soil_h2o_1500kpa_200cm.tif')
    solvar <- terra::rast(paste(geodir, '\\soils\\', soil.files, sep=''))
    names(solvar) <- c('Soil.pH.0cm', 'Soil.pH.mean', 'Soil.pH.200cm',
                       'Soil.h2o.33.0cm', 'Soil.h2o.33.mean', 'Soil.h2o.33.200cm',
                       'Soil.h2o.1500.0cm', 'Soil.h2o.1500.mean', 'Soil.h2o.1500.200cm')
    solcl <- data.frame(var=names(solvar), cluster='Soil 1')
  }

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
