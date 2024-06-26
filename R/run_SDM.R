# function(spname, #name of species. this is the only required parameter, all else are optional
#          data, merge, #data will allow user to supply own points. merge defines whether or not they should be merged with latest BIEN points
#          threshold, f.score, #threshold creates a thresholded version of the output. f.score allows user to tune this threshold to their use
#          envi, var.select, #envi allows user to supply own environmental data. var.select allows for turning off variable selection algorithm
#          bbox, output) #bbox allows user to set own extent, default is sized by extent of the data. output allows for saving outputs locally
# require(biomod2)
# require(plyr)
# require(stringr)
# require(terra)
#' @export

run_SDM <- function(spname, domain=world(), res){
  #### 1.0 Load Environmental data and species data ####
  #if(sources(domain)==sources(world())){
  domain <- pops.sdm::l48()
  domain <- pops.sdm::state(c('Oregon', 'California'))
  domain <- pops.sdm::county(state='Oregon', names=c('Coos', 'Curry', 'Douglas', 'Jackson', 'Josephine'))+
    pops.sdm::county(state='California', names=c('Del Norte', 'Siskiyou', 'Humboldt', 'Trinity', 'Shasta',
                                                 'Plumas', 'Mendocino', 'Tehama', 'Glenn', 'Butte', 'Lake',
                                                 'Colusa', 'Sutter', 'Yuba', 'Sierra', 'Nevada', 'Placer',
                                                 'Sonoma', 'Napa', 'Yolo', 'Sacramento', 'El Dorado', 'Marin',
                                                 'Amador', 'Alpine','Solano', 'Contra Costa', 'Calaveras',
                                                 'Tuolumne', 'San Francisco', 'San Mateo', 'San Joaquin',
                                                 'Alameda',  'Merced', 'Stanislaus', 'Santa Cruz', 'Mariposa',
                                                 'Santa Clara', 'Monterey', 'San Benito', 'San Luis Obispo',
                                                 'Santa Barbara', 'Ventura')); res <- 33

  domain <- pops.sdm::getEOR(); res <- 1000
  spname <- "Notholithocarpus densiflorus"

  myName <- stringr::str_replace(tolower(spname),' ', '_')
  dir <- getwd()

  #### 1.1 Gather Environmental Data ####
  #envi.vars <- pops.sdm::get_Envi(bio=T, lc=T, ptime=F, soil=T, pop=F, elev=T, res=res)
  envi <- pops.sdm::get_Envi(bio=T, elev=T, gdd=F, lc=T, pop=F, ptime=F, rnr=F, soil=F, tbase=5, res=res)
  base.r <- terra::crop(x=pops.sdm::rasterbase(res=res), y=domain, mask=T)
  base.r <- terra::subst(base.r, from=0, to=NA)

  if(res>250){
    envi.r <- terra::crop(x=envi[[1]], y=domain, mask=T)
    envi.r <- envi.r*base.r
    envi.cl <- pops.sdm::get_Clusters(envi.r)
    envi.vars <- list(envi.r, envi.cl$cluster); names(envi.vars) <- c('rast', 'clust')
  }
  if(res<=250){
    if(!file.exists(paste(dir, '/envi/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi[[1]]), '.tif', sep=''))){
      i <- 1
      for(i in i:terra::nlyr(envi[[1]])){print(i)
        envi.i <- terra::crop(x=envi[[1]][[i]], y=terra::ext(domain))
        if(i==1){envi.r <- envi.i}
        if(i>1){envi.r <- c(envi.r, envi.i)}
      }
      envi.r <- envi.r*base.r

      if(grepl('land', names(envi))){
      lvls.all <- data.frame(id=c(12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
                             landcover=as.factor(c('Ice', 'Dev_1', 'Dev_2', 'Dev_3', 'Dev_4', 'Barren', 'Decid', 'Everg',
                                                   'Mixed', 'Shrub', 'Grass', 'Pastr', 'Culti', 'WetWdy', 'WetHrb')))
      lcvar <- which(names(envi[[1]])=='landcover')
      envi[[1]][[lcvar]] <- terra::categories(envi[[1]][[lcvar]], value=lvls.all)
      }

      envi.cl <- pops.sdm::get_Clusters(envi.r)

      if(grepl('land', names(envi))){
        names(envi.tr$cluster)[which(is.na(names(envi.tr$cluster)))] <- 'landcover'
        envi.tr$cluster[which(names(envi.tr$cluster)=='landcover')] <- max(envi.tr$cluster)+1
      }

      envi.vars <- list(envi.r, envi.cl$cluster); names(envi.vars) <- c('rast', 'clust')
      terra::writeRaster(envi.r, paste(dir, '/envi/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))
    }
    if(file.exists(paste(dir, '/envi/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))){
      envi.r <- terra::rast(paste(dir, '/envi/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi[[1]]), '.tif', sep=''))
    }
  }

  #### 1.2 Gather Points ####
  pts.1 <- pops.sdm::get_pts.1(spname=spname, domain=domain)
  pts.r <- terra::rasterize(x=pts.1, y=base.r, fun='length', background=0)
  pts.r <- (pts.r*(base.r))>0
  pts.v <- terra::values(pts.r)
  pts.t <- which(pts.v==1)
  pts.f <- which(pts.v==0)
  pts.s <- sample(x=pts.f, size=length(pts.t))
  #myResp <- pts.v[c(pts.t, pts.s)]; myResp[myResp==0] <- NA
  myXY <- terra::xyFromCell(object=pts.r, cell=c(pts.t, pts.s))

  # if("landcover"%in%names(envi.vars[[1]])){
  #   myExtr <- terra::extract(envi.r$landcover, myXY)
  #   lvls.in <- unique(myExtr$landcover)
  #   lvls.all <- data.frame(id=c(12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
  #                          landcover=as.factor(c('Ice', 'Dev_1', 'Dev_2', 'Dev_3', 'Dev_4', 'Barren', 'Decid', 'Everg',
  #                                                'Mixed', 'Shrub', 'Grass', 'Pastr', 'Culti', 'WetWdy', 'WetHrb')))
  #   if(length(lvls.in)<nrow(lvls.all)){
  #     lvls.no <- lvls.all$id[which((lvls.all$id%in%lvls.in)==F)]
  #     envi.r$landcover <- terra::subst(envi.r$landcover, from=lvls.no, to=NA)
  #     lvls.all <- lvls.all[which(lvls.all$id!=lvls.no),]
  #   }
  #   envi.vars[[1]] <- terra::categories(envi.vars[[1]], value=lvls.all, layer=which(names(envi.r)=="landcover"))
  # }

  best.dir <- paste(dir, '/envi/envibest.', myName, '.', res, '.', names(envi), '.', terra::nrow(base.r), terra::ncol(base.r), '.tif', sep='')
  #### 2.0 Run variable selection function to choose the ideal variables ####
  if(!file.exists(best.dir)){
    envi.best <- pops.sdm::get_BestVars(envi=envi.vars$rast, pts=pts.r, clust=envi.vars$clust)
    terra::writeRaster(envi.best, best.dir)
  }
  if(file.exists(best.dir)){
    envi.best <- raster::stack(terra::rast(best.dir))
  }

  pts.2 <- pts.r!=is.na(pts.r)
  myResp <- terra::values(pts.2)
  myPA <- ifelse(myResp == 1, TRUE, FALSE)
  myPAtable <- data.frame(PA1 = myPA, PA2 = myPA, PA3 = myPA)
  pa.x <- vector()
  for (i in 1:ncol(myPAtable)){
    if(length(pts.t)<3333){npa <- length(pts.t)*3}
    if(length(pts.t)>3333&&length(pts.t)<=10000){npa <- 10000}
    if(length(pts.t)>10000){npa <- length(pts.t)}
    pa.s <- sample(which(myPAtable[, i] == FALSE), npa)
    myPAtable[pa.s, i] <- TRUE
    pa.x <- c(pa.x, pa.s)
  }
  pa.x <- unique(pa.x)
  myResp <- myResp[c(pts.t, pa.x)]
  myResp.PA <- ifelse(myResp == 1, 1, NA)
  myResp.XY <- terra::xyFromCell(object=pts.2, cell=c(pts.t, pa.x))
  myPAtable <- data.frame(myPAtable[c(pts.t, pa.x),])


  #### 3.0 Define options and parameters for modeling ####
  #mod.dir <- 'C:\Users\bjselige\Documents\pops.sdm\notholithocarpus.densiflorus'
  myExpl <- envi.best
  myOptions <- biomod2::BIOMOD_ModelingOptions()
  myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp.PA,
                                          expl.var = myExpl,
                                          resp.name = 'test',
                                          resp.xy = myResp.XY,
                                          PA.strategy = 'user.defined',
                                          PA.user.table = myPAtable)
  # Notes on algorithm choices; CTA is redundant with Random Forest, FDA and SRE have relatively low performance
  myAlgos <- c('GAM', 'GBM', 'GLM', 'RF', 'ANN', 'MARS', 'MAXENT')
  # Notes on evaluation methods: POD/SR/FR/BIAS is not useful, KAPPA, and TSS get similar results
  myEvals <- c('TSS') #myEvals <- c('ACCURACY', 'CSI', 'ETS', 'ROC', 'TSS')


  #### 4.0 Run and evaluate the models ####
  myModels <- biomod2::BIOMOD_Modeling(bm.format = myData,
                                       modeling.id = paste(myName,"AllModels",sep=""),
                                       models = myAlgos,
                                       bm.options = myOptions,
                                       # CV.strategy = 'random',
                                       CV.nb.rep = 3,
                                       CV.perc = .75,
                                       # CV.strategy='kfold',
                                       # CV.k=4,
                                       metric.eval = myEvals,
                                       var.import = 0,
                                       #SaveObj = T,
                                       CV.do.full.models = F) # use this option if you don't want a data split
  myEval <- biomod2::get_evaluations(myModels) # get all models evaluation


  #### 5.0 projection over the globe under current conditions ####
  myProj <- biomod2::BIOMOD_Projection(bm.mod = myModels,
                                       new.env = envi.best,
                                       # new.env = envi.df,
                                       # new.env.xy = envi.xy,
                                       proj.name = 'current',
                                       models.chosen = myModels@models.computed,
                                       #compress = 'xz',
                                       build.clamping.mask = F,
                                       binary.meth = NULL)
  myProj2 <- biomod2::get_predictions(myProj) # if you want to make custom plots, you can also get the projected map


  #### 6.0 Ensembling the models ####
  myEnsemble <- biomod2::BIOMOD_EnsembleModeling(bm.mod =  myModels,
                                                 models.chosen = myModels@models.computed, #myModels@models.computed[which(myModels@models.computed%in%mysub)],
                                                 em.by = 'all', #'all' combines all algos and PA runs into a single ensemble.
                                                 metric.eval = 'TSS', #only mean.weight and committee averaging use this argument, leaving it as 'TSS' is faster
                                                 metric.select = 'TSS',
                                                 em.algo = c('EMwmean', "EMcv", 'EMci', 'EMca'))
  myEvalEM <- biomod2::get_evaluations(myEnsemble) # get evaluation scores

  myProjEM <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myEnsemble,
                                                  bm.proj = myProj,
                                                  #selected.models = 'all',
                                                  #binary.meth = NULL,
                                                  #nb.cpu=parallel::detectCores()/2, #doesn't work on windows
                                                  #compress = 'xz',
                                                  metric.filter='TSS',
                                                  build.clamping.mask = F)

  ##### 5. Thresholding
  # myBinary <- plyr::llply(.data=c(1:length(myEvals)),
  #                         .fun=function(X){X.Eval <- myEvals[[X]]
  #                         X.Binary <- raster::stack(paste(tolower(stringr::str_replace(spname, ' ', '.')),
  #                                                         '\\proj_current\\proj_current_',
  #                                                         tolower(stringr::str_replace(spname, ' ', '.')),
  #                                                         '_', X.Eval, 'bin.grd', sep=''))
  #                         return(X.Binary)}); myBinary <- raster::stack(unlist(myBinary))
  #
  # myBinaryEM <- raster::stack(paste(tolower(stringr::str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
  #                                   tolower(stringr::str_replace(spname, ' ', '.')), '_ensemble_', myEvals[[1]], 'bin.grd', sep=''))

  p2 <- terra::unwrap(myProjEM@proj.out@val)[[1]]; p2.min <- min(terra::values(p2),na.rm=T); p2.max <- max(terra::values(p2),na.rm=T)
  trs <- seq(p2.min, p2.max, by=((p2.max-p2.min)/100))

  trs.metrics <- plyr::ldply(.data=c(1:length(trs)),
                             .fun=function(X){
                               p2.tr <- p2>trs[X]
                               x.zero <- sum(terra::extract(x=p2.tr, y=pts.1)==0, na.rm=T)/length(pts.1)
                               x.area <- sum(terra::values(p2.tr), na.rm=T)/length(na.omit(p2[]))
                               x.score <- 1 - x.zero - x.area
                               df.out <- data.frame('trs'=trs[X], 'zero'=x.zero, 'area'=x.area, 'score'=x.score, row.names=trs[X])
                               return(df.out)}, .progress = 'text')

  tr.best <- trs.metrics[which.max(trs.metrics$score),]
  p.tr <- p2[[1]]>tr.best$trs


  #### Finishing and outputs
  myProjEM2 <- terra::unwrap(myProjEM@proj.out@val)
  EM.mean <- myProjEM2$test_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo
  EM.cv <- myProjEM2$test_EMcvByTSS_mergedData_mergedRun_mergedAlgo
  EM.sd <- (EM.cv*EM.mean)/100 #EM.sd <- raster::calc(x=myProj@proj@val, fun=sd, na.rm=T)
  EM.ciInf <- myProjEM2$test_EMciInfByTSS_mergedData_mergedRun_mergedAlgo
  EM.ciSup <- myProjEM2$test_EMciSupByTSS_mergedData_mergedRun_mergedAlgo
  EM.ca <- myProjEM2$test_EMcaByTSS_mergedData_mergedRun_mergedAlgo
  p.out <- c(EM.mean, EM.cv, EM.sd, EM.ca, EM.ciInf, EM.ciSup)
  names(p.out) <- c('mean', 'cv', 'sd', 'ca', 'ciInf', 'ciSup')

  p3 <- c((p.out$mean), (p.tr)); names(p3) <- c('Raw', 'Binary')
  myEval2 <- list('All'=myEval, 'Ensemble'=myEvalEM, 'trs.metrics'=list('best'=tr.best$trs, 'full'=trs.metrics))
  outlist <- list('Data'=myData, 'Model'=myModels, 'Evaluation'=myEval2, 'Projections'=list('All'=myProj2, 'Ensemble'=p3))

  meta.df <- cbind(list(spname=spname, #extent='extent',
                        resolution=res,
                        n.points=length(pts.1),
                        n.psuedoabs=npa,
                        n.pa.reps=ncol(myPAtable),
                        n.cv.reps=length(unique(myEval$run)),
                        sources=c('GBIF', 'BIEN'),
                        best.vars=names(envi.best),
                        all.vars=names(envi.vars$rast),
                        algorithmns=myAlgos,
                        eval.metrics=myEvals,
                        #weights='weights',
                        threshold=tr.best$trs))
  write.csv(meta.df, paste(dir, '/', spname, 'meta.df.csv', sep=''))
  terra::writeRaster(p.out, paste(dir, '/', spname, '.output.tif', sep=''))
  return(outlist)
}


# splist <- c(
#   'Ailanthus altissima', #treeofheaven
#   'Buxus', #Boxwood
#   'Juglans nigra'# Black walmnut
#   # 'Lonicera hispidula', # Honey suckle
#   # 'Notholithocarpus densiflorus', # Tanoak
#   # 'Pseudotsuga menziesii',
#   # 'Tamarix chinensis',
#   # 'Tsuga canadensis',
#   # 'Tsuga caroliniana',
#   # 'Umbellularia californica' # Bay laurel
# )
# sp_ems <- llply(.data=c(1:length(splist)),
#                 .fun=function(x){EM_mod(splist[[x]], usa=T)},
#                 .progress='text'); names(sp_ems) <- splist
# pts <- BIEN::BIEN_occurrence_species(species=spname)
# pts <-  read.csv('C:\\Users\\bjselige\\Documents\\tree_of_heaven\\Ailanthus.BIEN.csv')
# pts <- gbif(genus=strsplit(x=spname, split=' ')[[1]][1], species=strsplit(x=spname, split=' ')[[1]][2],
#             geo=T, sp=T, removeZeros=T, download=T, ext=extent(biovars[[1]]))
#source('C:\\Users\\bjselige\\host_map\\get_pts.1.R')
#pts <- get_pts.1(spname)
# myBinary <-stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
#                        tolower(str_replace(spname, ' ', '.')), '_', myEvals, 'bin.grd', sep=''))
# names(myBinary) <- rep(paste(myEvals, substr(myModels@models.computed, 15, 27), sep='_'), length(myEvals))
# myBinary <- stack(list.files('C:\\Users\\bjselige\\host_map\\juglans.nigra\\proj_current\\individual_projections\\', 'bin.grd', full=T))
# names(myBinary) <- rep(paste(myEvals, substr(myModels@models.computed, 15, 27), sep='_'), length(myEvals))
# #### output plot
# borders <- usa
# rast <- myProjEM@proj@val
# rpts <- rasterToPoints(rast)
# rdf <- as.data.frame(rpts)
# ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) +
#   geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) +
#   scale_fill_continuous(type='viridis') +
#   theme_void() + theme(legend.position='none')
# png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\usa.',
#           gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''),
#     height=1080, width=2160); plot(ggsdm); dev.off()
