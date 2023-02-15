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

run_SDM <- function(spname, domain=world(), res=1){

  #### 1.0 Load Environmental data and species data ####
  #if(sources(domain)==sources(world())){
  envi.1k <- pops.sdm::get_Envi1k(bio=T, elev=T, soil=T) #envi.cv <- envi.1k$clust
  envi <- terra::crop(x=envi.1k$rast, y=terra::vect(domain), mask=T)#; envi.1k$rast <- envi
  pts.1 <- pops.sdm::get_pts.1(spname=spname, domain=domain)

  if(res>1){envi <- terra::aggregate(envi, fact=res, fun='mean')}
  if(res<1){envi <- terra::disagg(envi, fact=1/res, method='bilinear')}

  #### 1.1 Prep data ####
  pts.r <- terra::rasterize(x=pts.1, y=envi, fun='length', background=0)
  pts.r <- (pts.r*(envi[[1]]*0+1))>0
  pts.2 <- terra::as.points(pts.r)
  myName <- stringr::str_replace(tolower(spname),' ', '_')
  myResp <- pts.2$lyr1 # the presence/absences data for our species
  myResp[myResp==0] <- NA # setting 'true absences' to NA
  myRespXY <- terra::crds(pts.2) # the XY coordinates of species data


  #### 2.0 Run variable selection function to choose the ideal variables ####
  envi.best <- pops.sdm::get_BestVars(list('rast'=envi, 'clust'=envi.1k$clust), pts.2)


  #### 3.0 Define options and parameters for modeling ####
  #mod.dir <- 'C:\Users\bjselige\Documents\pops.sdm\notholithocarpus.densiflorus'
  myOptions <- biomod2::BIOMOD_ModelingOptions()
  myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                          expl.var = envi.best,
                                          resp.xy = myRespXY,
                                          resp.name = myName,
                                          #dir.name =
                                          PA.nb.rep = 1,
                                          PA.strategy = 'random',
                                          PA.nb.absences = sum(myResp, na.rm=T))
  # Notes on algorithm choices; CTA is redundant with Random Forest, FDA and SRE have relatively low performance
  myAlgos <- c('GLM', 'GAM', 'GBM', 'RF', 'ANN')#, 'MAXENT.Phillips')
  # Notes on evaluation methods: POD/SR/FR is not useful, KAPPA, ACCURACY, TSS, and ETS all get the same results
  myEvals <- c('CSI', 'ROC', 'TSS', 'ACCURACY', 'ETS', 'BIAS')


  #### 4.0 Run and evaluate the models ####
  myModels <- biomod2::BIOMOD_Modeling(bm.format = myData,
                                       modeling.id=paste(myName,"AllModels",sep=""),
                                       models = myAlgos,
                                       bm.options = myOptions,
                                       nb.rep = 1, #3, #number of runs
                                       data.split.perc = 80,
                                       metric.eval = myEvals,
                                       var.import = 1,
                                       #prevalence = 0.5,
                                       #save.output = T, # recommended to leave true
                                       do.full.models = F # use this option if you don't want a data split
  )
  myEval <- biomod2::get_evaluations(myModels) # get all models evaluation

  #### 5.0 projection over the globe under current conditions ####
  myProj <- biomod2::BIOMOD_Projection(bm.mod = myModels,
                                       new.env = envi.best,
                                       proj.name = 'current',
                                       #models.chosen = 'all',
                                       #models.chosen = c('GLM', 'GAM', 'GBM', 'RF', 'ANN'),
                                       #metric.binary = myEvals,
                                       #compress = 'xz',
                                       #build.clamping.mask = F,
                                       nb.cbu=32#,
                                       #output.format = '.grd',
                                       #do.stack=T
  )
  myProj2 <- biomod2::get_predictions(myProj) # if you want to make custom plots, you can also get the projected map


  #### 6.0 Ensembling the models ####
  myEnsemble <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myModels,
                                                 models.chosen = 'all',
                                                 em.by = 'PA_dataset+repet',
                                                 metric.eval = myEvals,
                                                 #metric.eval = NULL,
                                                 #eval.metric.quality.threshold = .5,
                                                 #models.eval.meth= myEvals,
                                                 prob.mean = T,
                                                 prob.cv = F, #don't use'
                                                 prob.ci = F, #prob.ci.alpha = 0.05,
                                                 prob.median = F,
                                                 committee.averaging = F,
                                                 prob.mean.weight = F, # leave on true
                                                 prob.mean.weight.decay = 'proportional')
  myEvalEM <- biomod2::get_evaluations(myEnsemble)[[1]] # get evaluation scores

  myProjEM <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myEnsemble,
                                                  bm.proj = myProj,
                                                  models.chosen = 'all',
                                                  metric.binary = myEvals,
                                                  compress = 'xz',
                                                  build.clamping.mask = F,
                                                  #output.format = '.grd'
  )

  myBinary <- plyr::llply(.data=c(1:length(myEvals)),
                          .fun=function(X){X.Eval <- myEvals[[X]]
                          X.Binary <- stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
                                                  tolower(str_replace(spname, ' ', '.')), '_', X.Eval, 'bin.grd', sep=''))
                          names(X.Binary) <- paste(X.Eval, substr(myProj@models.projected, nchar(spname)+1, max(nchar(myProj@models.projected))), sep='')
                          return(X.Binary)}); myBinary <- stack(unlist(myBinary))

  myBinaryEM <- raster::stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
                                    tolower(str_replace(spname, ' ', '.')), '_ensemble_', myEvals[[1]], 'bin.grd', sep=''))
  names(myBinaryEM) <- myEvals
  p.out <- myProjEM@proj@val[[1]]

  ##### 5. Thresholding
  p2 <- p.out
  trs <- seq(min(values(p2),na.rm=T), max(values(p2),na.rm=T), by=(max(values(p2),na.rm=T)-min(values(p2),na.rm=T))/100)
  trs.metrics <- plyr::ldply(.data=c(1:length(trs)),
                             .fun=function(X){
                               p2.tr <- rast(p2>trs[X])
                               x.zero <- sum(extract(x=p2.tr, y=pts)==0, na.rm=T)/length(pts)
                               x.area <- sum(values(p2.tr), na.rm=T)/length(na.omit(p2[]))
                               x.score <- 1 - x.zero - x.area
                               df.out <- data.frame('trs'=trs[X], 'zero'=x.zero, 'area'=x.area, 'score'=x.score, row.names=trs[X])
                               return(df.out)
                             }, .progress = 'text')

  tr.best <- trs.metrics[which.max(trs.metrics$score),]
  p.tr <- p.out>tr.best$trs
  p3 <- raster::stack(myProjEM@proj@val[[1]], p.tr); names(p3) <- c('Raw', 'Binary')
  p3 <- raster::stack(p3, myBinaryEM)
  myEval2 <- list(Eval=myEval, Ensemble=myEvalEM, BJS.metrics=list('best'=tr.best, 'full'=trs.metrics))
  outlist <- list(myData, myModels, myEval2, myProj2, p3)
  names(outlist) <- c('Data', 'Model', 'Evaluation', 'Projections', 'Ensemble')
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
