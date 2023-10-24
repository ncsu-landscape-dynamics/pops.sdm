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
  #domain <- pops.sdm::l48()
  domain <- pops.sdm::state(c('Oregon'))
  spname <- "Notholithocarpus densiflorus"; myName <- stringr::str_replace(tolower(spname),' ', '_')
  res <- 33
  dir <- getwd()


  #### 1.1 Gather Environmental Data ####
  envi.vars <- pops.sdm::get_Envi1k(bio=T, lc=T, ptime=F, soil=F, pop=F, res=res)
  base.r <- terra::crop(x=pops.sdm::rasterbase(res=res), y=domain, mask=T)

  if(!file.exists(paste(dir, '/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))){
    envi.r <- terra::crop(x=envi.vars$rast, y=domain, mask=T)
    envi.r <- envi.r*base.r
    terra::writeRaster(envi.r, paste(dir, '/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))
  }
  if(file.exists(paste(dir, '/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))){
    envi.r <- terra::rast(paste(dir, '/envi.', terra::nrow(base.r), terra::ncol(base.r), terra::nlyr(envi.vars$rast), '.tif', sep=''))
  }

  #### 1.2 Gather Points ####
  pts.1 <- pops.sdm::get_pts.1(spname=spname, domain=domain)
  # myResp <- rep(1, length(pts.1))#pts.2 <- sp::SpatialPoints(terra::crds(pts.1))
  # myXY <- terra::crds(pts.1)
  # PA.df <-

  #if(!file.exists()){
  pts.r <- terra::rasterize(x=pts.1, y=base.r, fun='length', background=0)
  pts.r <- (pts.r*(base.r))>0
  pts.v <- terra::values(pts.r)
  pts.t <- which(pts.v==1)
  pts.f <- which(pts.v==0)
  pts.s <- sample(x=pts.f, size=length(pts.t))
  myResp <- pts.v[c(pts.t, pts.s)]; myResp[myResp==0] <- NA
  myXY <- terra::xyFromCell(object=pts.r, cell=c(pts.t, pts.s))

  # #pts.2 <- terra::as.points(pts.r); names(pts.2) <- 'lyr1'
  # #myResp <- pts.2[[1]]; # the presence/absences data for our species
  # pts.t <- which(pts.2$lyr1==1)
  # pts.f <- which(pts.2$lyr1==0)
  # pts.r <- sample(x=pts.f, size=length(pts.t))
  # pts.3 <- pts.2[c(pts.t, pts.r)]
  # myResp <- pts.3$lyr1
  # myResp[myResp==0] <- NA # setting 'true absences' to NA
  # myXY <- terra::crds(pts.3)#myXY <- terra::crds(pts.2) # the XY coordinates of species data
  #
  # PA.df <- as.data.frame(myResp); PA.df[is.na(PA.df)] <- FALSE
  # PA.fact <- sum(PA.df==F)/sum(myResp, na.rm=T)
  # PA.df$myResp[which(PA.df==F)[round((1:sum(myResp, na.rm=T))*PA.fact)]] <- TRUE
  # PA.df$myResp <- as.logical(PA.df$myResp)
  #}

  #### 2.0 Run variable selection function to choose the ideal variables ####
  envi.best <- pops.sdm::get_BestVars(list('rast'=envi.r, 'clust'=envi.vars$clust), pts.2)
  envi.best <- raster::stack(envi.best)

  #### 3.0 Define options and parameters for modeling ####
  #mod.dir <- 'C:\Users\bjselige\Documents\pops.sdm\notholithocarpus.densiflorus'
  myOptions <- biomod2::BIOMOD_ModelingOptions()
  myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                          expl.var = envi.best,
                                          resp.xy = myXY,
                                          resp.name = myName,
                                          PA.nb.rep = 1,
                                          PA.strategy = 'random', #, 'user.defined', #
                                          PA.nb.absences = sum(myResp, na.rm=T)
                                          #PA.table = PA.df
                                          )
  # Notes on algorithm choices; CTA is redundant with Random Forest, FDA and SRE have relatively low performance
  myAlgos <- c('GAM', 'GBM', 'GLM', 'RF', 'ANN', 'MARS', 'MAXENT.Phillips')
  # Notes on evaluation methods: POD/SR/FR/BIAS is not useful, KAPPA, and TSS get similar results
  myEvals <- c('ACCURACY', 'CSI', 'ETS', 'ROC', 'TSS')


  #### 4.0 Run and evaluate the models ####
  myModels <- biomod2::BIOMOD_Modeling(data = myData,
                                       modeling.id = paste(myName,"AllModels",sep=""),
                                       models = 'GLM',
                                       models.options = myOptions,
                                       NbRunEval = 1, #3, #number of runs
                                       DataSplit = 80,
                                       models.eval.meth = myEvals,
                                       VarImport = 1,
                                       SaveObj = T,
                                       do.full.models = F) # use this option if you don't want a data split
  myEval <- biomod2::get_evaluations(myModels) # get all models evaluation


  #### 5.0 projection over the globe under current conditions ####
  myProj <- biomod2::BIOMOD_Projection(modeling.output = myModels,
                                       new.env = envi.best,
                                       proj.name = 'current',
                                       #compress = 'xz',
                                       #build.clamping.mask = F,
                                       #output.format = '.grd',
                                       #do.stack=T
                                       #nb.cbu=32,
                                       binary.meth = myEvals)
  myProj2 <- biomod2::get_predictions(myProj) # if you want to make custom plots, you can also get the projected map

  vals1 <- na.omit(values(envi.best))
  vals2 <- na.omit(values(myProj@proj@val))

  myname2 <- stringr::str_replace(myName, '_', '.')
  mysub <- c(paste(myname2, '_PA1_RUN1_GAM', sep=''),
             paste(myname2,'_PA1_RUN1_GLM', sep=''),
             paste(myname2,'_PA1_RUN1_RF', sep=''),
             paste(myname2,'_PA1_RUN1_ANN', sep=''),
             paste(myname2,'_PA1_RUN1_MARS', sep=''))
  #### 6.0 Ensembling the models ####
  myEnsemble <- biomod2::BIOMOD_EnsembleModeling(modeling.output =  myModels,
                                                 chosen.models = myModels@models.computed[which(myModels@models.computed%in%mysub)],
                                                 em.by = 'all', #'all' combines all algos and PA runs into a single ensemble.
                                                 eval.metric = 'TSS', #only mean.weight and committee averaging use this argument, leaving it as 'TSS' is faster
                                                 models.eval.meth = myEvals,
                                                 prob.mean = T,
                                                 prob.cv = T,
                                                 prob.ci = T, prob.ci.alpha = 0.05,
                                                 #prob.median = F, prob.mean.weight = T,
                                                 #prob.mean.weight.decay = 'proportional',
                                                 committee.averaging = F)

  myEnsembleCA <- biomod2::BIOMOD_EnsembleModeling(modeling.output = myModels,
                                                   chosen.models = myModels@models.computed[which(myModels@models.computed%in%mysub)],
                                                   em.by = 'all', #'all' combines all algos and PA runs into a single ensemble.
                                                   eval.metric = 'all',#   eval.metric = 'TSS',
                                                   models.eval.meth = myEvals,
                                                   prob.mean = F, prob.cv = F, prob.ci = F, #prob.ci.alpha = 0.05,
                                                   committee.averaging = T)
  myEvalEM <- biomod2::get_evaluations(myEnsemble) # get evaluation scores
  myEvalCA <- biomod2::get_evaluations(myEnsembleCA)

  myProjEM <- biomod2::BIOMOD_EnsembleForecasting(EM.output = myEnsemble,
                                                  projection.output = myProj,
                                                  selected.models = 'all',
                                                  binary.meth = myEvals,
                                                  compress = 'xz',
                                                  build.clamping.mask = F)
  myProjCA <- biomod2::BIOMOD_EnsembleForecasting(EM.output = myEnsembleCA,
                                                  projection.output = myProj,
                                                  selected.models = 'all',
                                                  binary.meth = myEvals,
                                                  compress = 'xz',
                                                  build.clamping.mask = F)
  EM.mean <- raster::mean(myProjEM@proj@val[[which(grepl('EMmean', names(myProjEM@proj@val)))]])
  EM.cv <- raster::mean(myProjEM@proj@val[[which(grepl('EMcv', names(myProjEM@proj@val)))]])
  EM.sd <- (EM.cv*EM.mean)/100 #EM.sd <- raster::calc(x=myProj@proj@val, fun=sd, na.rm=T)
  EM.ciInf <- raster::mean(myProjEM@proj@val[[which(grepl('EMciInf', names(myProjEM@proj@val)))]])
  EM.ciSup <- raster::mean(myProjEM@proj@val[[which(grepl('EMciSup', names(myProjEM@proj@val)))]])
  EM.ca <- raster::mean(myProjCA@proj@val)
  p.out <- raster::stack(EM.mean, EM.cv, EM.sd, EM.ca, EM.ciInf, EM.ciSup)
  names(p.out) <- c('mean', 'cv', 'sd', 'ca', 'ciInf', 'ciSup')
  #p.out <- mean(terra::unwrap(myProjEM@proj.out@val))


  myBinary <- plyr::llply(.data=c(1:length(myEvals)),
                          .fun=function(X){X.Eval <- myEvals[[X]]
                          X.Binary <- raster::stack(paste(tolower(stringr::str_replace(spname, ' ', '.')),
                                                          '\\proj_current\\proj_current_',
                                                          tolower(stringr::str_replace(spname, ' ', '.')),
                                                          '_', X.Eval, 'bin.grd', sep=''))
                          #names(X.Binary) <- paste(X.Eval, substr(myProj@models.projected, nchar(spname)+1, max(nchar(myProj@models.projected))), sep='')
                          return(X.Binary)}); myBinary <- raster::stack(unlist(myBinary))

  myBinaryEM <- raster::stack(paste(tolower(stringr::str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
                                    tolower(stringr::str_replace(spname, ' ', '.')), '_ensemble_', myEvals[[1]], 'bin.grd', sep=''))
  #names(myBinaryEM) <- myEvals

  ##### 5. Thresholding
  p2 <- p.out$mean
  p2.min <- min(terra::values(p2),na.rm=T); p2.max <- max(terra::values(p2),na.rm=T)
  trs <- seq(p2.min, p2.max, by=((p2.max-p2.min)/100))

  trs.metrics <- plyr::ldply(.data=c(1:length(trs)),
                             .fun=function(X){
                               p2.tr <- p2>trs[X]
                               x.zero <- sum(raster::extract(x=p2.tr, y=as(pts.1, 'Spatial'))==0, na.rm=T)/length(pts.1)
                               x.area <- sum(terra::values(p2.tr), na.rm=T)/length(na.omit(p2[]))
                               x.score <- 1 - x.zero - x.area
                               df.out <- data.frame('trs'=trs[X], 'zero'=x.zero, 'area'=x.area, 'score'=x.score, row.names=trs[X])
                               return(df.out)}, .progress = 'text')

  tr.best <- trs.metrics[which.max(trs.metrics$score),]
  p.tr <- p.out$mean>tr.best$trs
  p3 <- raster::stack((p.out$mean), (p.tr)); names(p3) <- c('Raw', 'Binary')#p3 <- raster::stack(p3, myBinaryEM)
  myEval2 <- list('All'=myEval, 'Ensemble'=myEvalEM, 'trs.metrics'=list('best'=tr.best$trs, 'full'=trs.metrics))
  outlist <- list('Data'=myData, 'Model'=myModels, 'Evaluation'=myEval2, 'Projections'=list('All'=myProj2, 'Ensemble'=p3))
  return(outlist)

  meta.df <- data.frame(spname=spname,
                        #extent='extent',
                        resolution=res,
                        n.points=length(pts.1),
                        sources=c('GBIF', 'BIEN'),
                        best.vars=names(envi.best),
                        all.vars=names(envi.vars),
                        algorithmns=myAlgos,
                        tss=myEvals,
                        weights='weights',
                        threshold=tr.best$trs)

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
