# require(geodata)
# require(terra)
# require(parallel)
#' @export

get_Envi1k <- function(bio=F, elev=F, gdd=F, lc=F, pop=F, ptime=F, rnr=F, soil=F, tbase=5, res=1000){

  geodir <- 'Q:\\Shared drives\\Data\\Raster\\'
  base.rast <- pops.sdm::rasterbase(res=res)

  if(bio==T){
    if(res>=1000){biovar <- geodata::worldclim_global(var='bio', res=.5, path=paste(geodir, 'Global\\',sep=''))}
    if(res<1000){
      if(!file.exists(paste(geodir, 'USA\\bioclim\\500m\\BioClimComposite_1971_2000_15s.tif', sep=''))){
        bio <- terra::rast('Q:\\Shared drives\\Data\\Original\\BioClimComposite_1971_2000_400m.tif')
        bio <- bio[[c(1:4, 6:20)]]
        bio.p <- terra::project(bio, base.rast, threads=T)
        terra::writeRaster(bio.m, filename=paste(geodir, 'USA\\bioclim\\500m\\BioClimComposite_1971_2000_15s.tif', sep=''))
      }
      if(file.exists(paste(geodir, 'USA\\bioclim\\500m\\BioClimComposite_1971_2000_15s.tif', sep=''))){
        biovar <- terra::rast(paste(geodir, 'USA\\bioclim\\500m\\BioClimComposite_1971_2000_15s.tif', sep=''))
      }
      if(res<500){biodir <- paste(geodir, 'USA\\bioclim\\', res, 'm\\', sep='')
      if(!dir.exists(biodir)){dir.create(biodir)}
      if(dir.exists(biodir)){
        if(length(list.files(biodir))<19){i <- 1
        for(i in i:terra::nlyr(biovar)){print(i)
          i.biovar <- terra::project(biovar[[i]], base.rast, threads=T)
          terra::writeRaster(i.biovar, filename=paste(biodir, 'bioclim.', i, '_', res, 'm.tif', sep=''))
        }
        }
        if(length(list.files(biodir))==19){
          bio.list <- list(); i <- 1
          for(i in i:19){bio.list[[i]] <- terra::rast(paste(biodir, 'bioclim.', i, '_', res, 'm.tif', sep=''))}
          biovar <- terra::rast(bio.list)
        }
      }
      }
      if(res>500){biovar <- terra::aggregate(biovar, fact=(res/500), fun='mean')} #defunct, no resolutions planned between 500 and 1000
    }
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
    if(res>=1000){
      elevvar <- geodata::elevation_global(res=.5, path=paste(geodir, 'Global\\',sep=''))
      names(elevvar) <- 'Elevation'
    }
    #   if(res<1000){}
    elevcl <- data.frame(var=names(elevvar), cluster='Elevation 1')
  }

  if(gdd==T){
    if(exists('tbase')==F){tbase <- 5}
    getGDD <- function(tbase){
      gdpath <- paste(paste(geodir, 'Global\\',sep=''), 'gdd.base', tbase, '.tif', sep='')
      if(!file.exists(gdpath)){print('Calculating GDD')
        tavg <- geodata::worldclim_global(var='tavg', res=.5, path=paste(geodir, 'Global\\',sep=''))
        g1 <- terra::app(x=tavg, tbase=tbase, fun=function(x, tbase){
          days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
          xbase <- x-tbase; xbase[which(xbase<0)] <- 0; xbase <- sum(xbase*days)
          return(xbase)}, cores=parallel::detectCores()/2); closeAllConnections()
        terra::writeRaster(g1, filename=gdpath)
      }
      if(file.exists(gdpath)){g1 <- terra::rast(gdpath)}

      if(res<1000){
        gdpath <- paste(geodir, 'USA\\gdd_base',tbase, '_', res, 'm.tif', sep='')
        if(!file.exists(gdpath)){
          g1 <- terra::project(g1, base.rast, threads=T)
          terra::writeRaster(g1, filename=gdpath)
        }
        if(file.exists(gdpath)){g1 <- terra::rast(gdpath)}
      }
      names(g1) <- paste('GDD ', tbase, sep=''); return(g1)
    }
    gddvar <- getGDD(tbase)
    gddcl <- data.frame(var=names(gddvar), cluster='Temp 1')
  }

  if(lc==T){
    if(res>=1000){
      lcpath <- paste(geodir, 'Global\\landcover\\', sep='')
      if(dir.exists(lcpath)){
        built <- terra::rast(paste(lcpath, 'built.tif', sep=''))
        cropl <- terra::rast(paste(lcpath, 'cropl.tif', sep=''))
        grass <- terra::rast(paste(lcpath, 'grass.tif', sep=''))
        shrub <- terra::rast(paste(lcpath, 'shrub.tif', sep=''))
        trees <- terra::rast(paste(lcpath, 'trees.tif', sep=''))
        wetld <- terra::rast(paste(lcpath, 'wetld.tif', sep=''))
      }
      if(!dir.exists(lcpath)){
        globe <- geodata::worldclim_global(var='bio', res=.5, path=geodir); globext <- terra::ext(globe)
        built <- terra::extend(x=geodata::landcover(var='built', path=geodir), y=globext)
        cropl <- terra::extend(x=geodata::landcover(var="cropland", path=geodir), y=globext)
        grass <- terra::extend(x=geodata::landcover(var='grassland', path=geodir), y=globext)
        shrub <- terra::extend(x=geodata::landcover(var='shrubs', path=geodir), y=globext)
        trees <- terra::extend(x=geodata::landcover(var='trees', path=geodir), y=globext)
        wetld <- terra::extend(x=geodata::landcover(var='wetland', path=geodir), y=globext)
        terra::writeRaster(built, filename=paste(lcpath, 'built.tif', sep=''))
        terra::writeRaster(cropl, filename=paste(lcpath, 'cropl.tif', sep=''))
        terra::writeRaster(grass, filename=paste(lcpath, 'grass.tif', sep=''))
        terra::writeRaster(shrub, filename=paste(lcpath, 'shrub.tif', sep=''))
        terra::writeRaster(trees, filename=paste(lcpath, 'trees.tif', sep=''))
        terra::writeRaster(wetld, filename=paste(lcpath, 'wetld.tif', sep=''))
        lcvar <- c(built, cropl, grass, shrub, trees, wetld)
      }
    }
    if(res<1000){
      lcpath <- paste(geodir, 'USA\\landcover\\', sep='')
      built <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_built_21_22_23_24.tif', sep=''))
      decid <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_decid_41_43.tif', sep=''))
      everg <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_everg_42_43.tif', sep=''))
      trees <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_trees_41_42_43.tif', sep=''))
      shrub <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_shrub_52.tif', sep=''))
      grass <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_grass_71.tif', sep=''))
      pastr <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_pastr_81.tif', sep=''))
      cropl <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_cropl_82.tif', sep=''))
      culti <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_culti_81_82.tif', sep=''))
      wetld <- terra::rast(paste(lcpath, '\\nlcd_2019_1s_wetld_90_95.tif', sep=''))

      lcvar <- c(built, decid, everg, trees, shrub, grass, pastr, cropl, culti, wetld)
      names(lcvar) <- c('built', 'decid', 'everg', 'trees', 'shrub', 'grass', 'pastr', 'cropl', 'culti', 'wetld')

      # if(!file.exists(paste(lcpath, 'nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))){
      #   lc30 <- terra::rast(paste(lcpath, 'nlcd_2019_land_cover_l48_20210604.tif', sep=''))
      #   lcbase <- rasterbase(res=33)
      #   lc1s <- terra::project(lc30, lcbase, method='near', threads=T)
      #   # lc30.p <- terra::project(lc30, "epsg:4326", method='near', threads=T)
      #   # lc30.c <- terra::crop(lc30.p, terra::ext(pops.sdm::l48()))
      #   # lc30.m <- terra::crop(lc30.m, y=pops.sdm::l48(), mask=T)
      #   terra::writeRaster(lc1s, paste(lcpath, 'nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))
      # }
      # if(file.exists(paste(geodir, 'USA\\landcover\\', 'nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))){
      #   lcvar <- terra::rast(paste(lcpath, 'nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))
      # }
      if(res>33){
        if(!file.exists(paste(lcpath, 'ncld_2019_', res, 'm.tif', sep=''))){
          if(res==250){
            lcagg <- terra::aggregate(lcvar, fact=7, fun='modal')
            lcagg <- terra::project(lcagg, base.rast, method='near', threads=T)
          }
          if(res!=250){
            lcagg <- terra::aggregate(lcvar, fact=(res/(100/3)), fun='modal')
          }
          terra::writeRaster(lcagg, paste(lcpath, 'ncld_2019_', res, 'm.tif', sep=''))
        }
        if(file.exists(paste(lcpath, 'ncld_2019_', res, 'm.tif', sep=''))){
          lcvar <- terra::rast(paste(lcpath, 'ncld_2019_', res, 'm.tif', sep=''))
        }
      }
    }
    lccl <- data.frame(var=names(lcvar), cluster=NA)
    lccl$cluster[which(lccl$var%in%c('built'))] <- 'Development'
    lccl$cluster[which(lccl$var%in%c('decid', 'everg', 'trees', 'shrub', 'grass', 'wetld'))] <- 'Vegetation'
    lccl$cluster[which(lccl$var%in%c('pastr', 'cropl', 'culti'))] <- 'Agriculture'
  }

  if(pop==T){
    if(res>=1000){popvar <- geodata::population(year='2020', res=.5, path=paste(geodir, 'Global\\',sep=''))}
    if(res<1000){
      if(!file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_3ss_V1_0.tif', sep=''))){
        popvar <- terra::rast(paste(geodir, 'Global\\pop\\GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif', sep=''))
        popvar.c <- terra::crop(popvar, terra::ext(pops.sdm::l48()))
        popvar.m <- terra::crop(x=popvar.c, y=pops.sdm::l48(), mask=T)
        terra::writeRaster(popvar.m, paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_3ss_V1_0.tif', sep=''))
      }
      if(file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_3ss_V1_0.tif', sep=''))){
        popvar <- terra::rast(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_3ss_V1_0.tif', sep=''))
      }
      if(res>100){
        if(!file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))){
          popagg <- terra::project(popvar, base.rast, threads=T) #popagg <- terra::aggregate(popvar, fact=(res/100), fun='mean')
          terra::writeRaster(popagg, paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))
        }
        if(file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))){
          popvar <- terra::rast(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))
        }
      }
      if(res<100){
        if(!file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))){
          popdis <- terra::project(popvar, base.rast, threads=T)
          terra::writeRaster(popdis, paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))
        }
        if(file.exists(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))){
          popvar <- terra::rast(paste(geodir, 'USA\\GHS_POP_E2020_USA_R2023A_4326_', res, 'm.tif', sep=''))
        }
      }
    }
    names(popvar) <- 'Population'; popcl <- data.frame(var=names(popvar), cluster='Development')
  }

  if(ptime==T){
    getPrecipTiming <- function(){
      if(res==1000){ptpath <- paste(geodir, 'Global\\precip_timing.tif', sep='')}
      if(res<1000){ptpath <- paste(geodir, 'USA\\precip_timing_', res, 'm.tif', sep='')}
      if(file.exists(ptpath)){prect <- terra::rast(ptpath)}
      if(!file.exists(ptpath)){print('Calculating Precip Timing (DJF-JJA)')
        prec <- geodata::worldclim_global(var='prec', res=.5, path=paste(geodir, 'Global\\',sep=''))
        if(res<1000){
          prec <- terra::project(prec, base.rast, threads=T)
        }
        prect <- terra::app(x=prec, cores=parallel::detectCores()/2,
                            fun=function(x){return(sum(x[[12]], x[[1]], x[[2]])-sum(x[[6]], x[[7]], x[[8]]))})
        terra::writeRaster(prect, filename=ptpath); closeAllConnections()
      }
      names(prect) <- 'Precip.Timing'; return(prect)
    }
    timevar <- getPrecipTiming()
    timecl <- data.frame(var=names(timevar), cluster='Precip 2')
  }

  if(rnr==T){
    rnrvar <- terra::rast(list.files(paste(geodir, 'Global\\',sep=''), '.dist.wrld', full.names = T))
    names(rnrvar) <- c('Road.Dist', 'Rail.Dist')
    rnrcl <- data.frame(var=names(rnrvar), cluster='Roads/Rails')
  }

  if(soil==T){
    if(res<1000){
      if(any(grepl('.mean', list.files(paste(geodir, 'USA\\soils\\250m\\', sep=''))))==F){
        if(any(grepl('.mean', list.files(paste(geodir, 'Global\\soils\\250m\\', sep=''))))==F){
          var <- c('ph.h2o', '1500kPa', '33kPa'); i <- 1
          for(i in i:length(type)){
            sl <- lapply(list.files(path=paste(geodir, 'Global\\soils\\250m\\', sep=''), pattern=var[i], full.names=T), terra::rast)
            sl.mean <- terra::mean(sl[[1]], sl[[2]], sl[[3]], sl[[4]], sl[[5]], sl[[6]])
            writeRaster(sl.mean, paste('C:\\Users\\bjselige\\Desktop\\Soil_OpenLandMap\\', var[i], '.mean.tif', sep=''))
          }
        }
        if(any(grepl('.mean', list.files(paste(geodir, 'Global\\soils\\250m\\', sep=''))))){
          file <- c('ph.mean', '1500kPa.mean', '33kPa.mean'); i <- 1
          for(i in i:length(file)){
            sl.mean <- terra::rast(paste(geodir, 'Global\\soils\\250m\\', file[i], '.tif', sep=''))
            sl.crop <- terra::crop(sl.mean, terra::ext(pops.sdm::l48()))
            sl.mask <- terra::crop(sl.crop, y=pops.sdm::l48(), mask=T)
            terra::writeRaster(sl.mask, paste(geodir, 'USA\\soils\\250m\\', file[i], '.tif', sep=''))
          }
        }
      }
      if(any(grepl('.mean', list.files(paste(geodir, 'USA\\soils\\250m\\', sep=''))))){
        solvar <- terra::rast(list.files(paste(geodir, 'USA\\soils\\250m\\', sep=''), pattern='.mean.', full.names=T))
        names(solvar) <- c('soil.1500kPa.mean', 'soil.33kPa.mean', 'soil.ph.mean')
      }

      if(res<250){
        soldir <- paste(geodir, 'USA\\soils\\', res, 'm\\', sep='')
        if(!dir.exists(soldir)){dir.create(soldir)}
        if(dir.exists(soldir)){
          if(length(list.files(soldir))<3){i <- 1
          for(i in i:terra::nlyr(solvar)){print(i)
            i.sol <- terra::project(solvar[[i]], base.rast, threads=T)
            terra::writeRaster(i.sol, filename=paste(soldir, names(solvar)[[i]], '_', res, 'm.tif', sep=''))
          }
          }
          if(!length(list.files(soldir))<3){solvar <- terra::rast(list.files(soldir, full.names=T))}
        }
      }

      if(res>250){
        soldir <- paste(geodir, 'USA\\soils\\', res, 'm\\', sep='')
        if(!dir.exists(soldir)){dir.create(soldir); i <- 1
        for(i in i:terra::nlyr(solvar)){
          i.solagg <- terra::aggregate(solvar[[i]], fact=(res/250), fun='mean')
          terra::writeRaster(i.solagg, filename=paste(soldir, names(solvar)[[i]], '_', res, 'm.tif', sep=''))
        }
        }
        if(dir.exists(soldir)){solvar <- terra::rast(list.files(soldir, 'soil.', full.names=T))}
      }

      names(solvar) <- c('soil.1500kPa.mean', 'soil.33kPa.mean', 'soil.ph.mean')

    }
    if(res>=1000){
      soil.files <- c('Soil_pH_0cm.tif', 'Soil_pH_mean.tif', 'Soil_pH_200cm.tif',
                      'Soil_h2o_33kpa_0cm.tif', 'Soil_h2o_33kpa_mean.tif',
                      'Soil_h2o_33kpa_200cm.tif', 'Soil_h2o_1500kpa_0cm.tif',
                      'Soil_h2o_1500kpa_mean.tif', 'Soil_h2o_1500kpa_200cm.tif')
      solvar <- terra::rast(paste(geodir, 'Global\\soils\\1km\\', soil.files, sep=''))
      names(solvar) <- c('Soil.pH.0cm', 'Soil.pH.mean', 'Soil.pH.200cm',
                         'Soil.h2o.33.0cm', 'Soil.h2o.33.mean', 'Soil.h2o.33.200cm',
                         'Soil.h2o.1500.0cm', 'Soil.h2o.1500.mean', 'Soil.h2o.1500.200cm')
    }
    solcl <- data.frame(var=names(solvar), cluster='Soil')
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
