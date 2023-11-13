# require(biomod2)
# require(plyr)
# require(raster)
# require(terra)
#' @export

get_BestVars <- function(envi, pts, clust){

  #if(class(envi.r)[1]=="SpatRaster"){envi2 <- raster::stack(envi$rast)}
  envi2 <- envi; envi.cv <- clust
  pts.v <- terra::values(pts)
  pts.t <- which(pts.v==1)
  pts.f <- which(pts.v==0)
  pts.s <- sample(x=pts.f, size=length(pts.t))
  myResp <- pts.v[pts.t]#; myResp[myResp==0] <- NA
  #myXY <- terra::xyFromCell(object=pts.r, cell=c(pts.t, pts.s))
  myXY <- terra::xyFromCell(object=pts.r, cell=pts.t)
  myPA <- terra::xyFromCell(object=pts.r, cell=pts.s)
  nreps <- 1 #unnecessary since running full models yields same results, no extra reps needed
  #p.max <- 1000000-length(pts.t) #pts.r <- sample(x=pts.f, size=pmin(length(pts.f), p.max))

  PA.df <- as.data.frame(myPA)

  # PA.df <- as.data.frame(myResp); PA.df[is.na(PA.df)] <- FALSE
  # PA.fact <- sum(PA.df==F)/sum(myResp, na.rm=T)
  # PA.df$myResp[which(PA.df==F)[round((1:sum(myResp, na.rm=T))*PA.fact)]] <- TRUE
  # PA.df$myResp <- as.logical(PA.df$myResp)


  # Notes on algorithm selection. GBM and ANN do not work with all data used/no split.
  # SRE, CTA, FDA, GLM, GAM, RF, and MARS are comparable in speed (SRE fastest, MARS slowest by about 25%, rest in respective order)
  # Maxent is approximately 2.5 to 3 times slower than the previous 7 algorithms mentioned
  # RF yields extremely high accuracy scores, may bias results, exchanged for CTA which is similar
  # Running all algorithms at once is approximately 25% faster than running them separately.
  #algos <- list('SRE', 'CTA', 'FDA', 'GLM', 'GAM', 'MARS', 'MAXENT.Phillips', c('SRE', 'CTA', 'FDA', 'GLM', 'GAM', 'MARS', 'MAXENT.Phillips'))
  #names(algos) <- c('SRE', 'CTA', 'FDA', 'GLM', 'GAM',  'MARS', 'MAXENT.Phillips', 'ALL')
  algos <- list(c('FDA', 'GLM', 'GAM', 'MARS', 'MAXENT.Phillips')); names(algos) <- c('ALL')
  evals <- c('ACCURACY', 'CSI', 'ETS', 'ROC', 'TSS') #Notes on evals; Kappa similar to tss, bias/far/sr/pod not very useful,

  i <- 1; t.list <- data.frame('algos'=names(algos), 'time'=rep(NA, length(algos)))
  k <- 1; k.list <- list()
  j <- 1; j.list <- list(); j.data <- data.frame()

  for(i in 1:length(algos)){i.t1 <- Sys.time()
    i.algo <- algos[[i]]; n.algo <- names(algos)[i];
    for(j in 1:nreps){#reps of same algo
      for(k in 1:max(envi.cv)){#reps through clusters
        if(k==1){k.names <- names(envi2); k.stack <-NULL}

        k.test <- plyr::ldply(.data=k.names,
                              .fun=function(X){

                                myExtr <- terra::extract(c(envi2[[X]], k.stack), myXY)

                                if("nlcd_2019_land_cover_l48_20210604"%in%colnames(myExtr)){
                                  myExtr$nlcd_2019_land_cover_l48_20210604 <- as.factor(myExtr$nlcd_2019_land_cover_l48_20210604)
                                }



                                myOptions <- biomod2::BIOMOD_ModelingOptions('GLM'=list(test='none'))
                                myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                                        expl.var = myExtr,
                                                                        resp.xy = myXY,
                                                                        resp.name = 'test',
                                                                        PA.nb.rep = 1,
                                                                        PA.strategy = 'user.defined',
                                                                        PA.user.table = PA.df)

                                # 3. Computing the models
                                myModels <- biomod2::BIOMOD_Modeling(bm.format = myData, #bm.format = myData,
                                                                     bm.options = myOptions, #bm.options = myOptions,
                                                                     models = i.algo,
                                                                     CV.nb.rep = 1, #nb.rep = 1, #number of runs
                                                                     CV.strategy = 'user.defined',
                                                                     CV.user.table = PA.df,
                                                                     CV.perc = 1,
                                                                     metric.eval = evals, # metric.eval = c('ROC', 'TSS'),
                                                                     # save.output = T, # recommended to leave true
                                                                     # scale.models = F, #experimental don't use
                                                                     CV.do.full.models = F,
                                                                     weights=NULL)

                                # Evaluation
                                options(digits = 22)
                                myEval <- biomod2::get_evaluations(myModels)
                                eval.v <- c(myEval[evals,"Testing.data",,,])
                                return(mean(eval.v))
                              })
        row.names(k.test) <- k.names

        k.out <- data.frame(cbind(envi.cv[row.names(k.test)], k.test$V1))
        colnames(k.out) <- c('clust', 'score')
        k.out <- k.out[order(k.out$clust),]
        k.out[,3] <- as.integer(rank(-k.out$score))
        k.var <- row.names(k.out)[which.min(k.out$V3)]
        k.score <- k.out$score[which.max(k.out$score)]
        k.clust <- k.out[k.var, 'clust']
        k.envi <- envi2[[k.var]]
        k.stack <- c(k.envi, k.stack)
        k.list[[k]] <- list('score'=k.score, 'stack'=k.stack)

        # this reduces the number of variables in each cluster to the top 2 (redux var) most important
        if(k==1){redux <- 2; k.redux <- data.frame()
        for(i in 1:length(unique(k.out$clust))){
          i.clust <- k.out[k.out$clust==i,]
          i.redux <- i.clust[order(i.clust$score, decreasing=T),][1:redux,]
          k.redux <- rbind(k.redux, i.redux)
        }
        k.out <- k.redux
        }
        k.names <-  row.names(k.out)[which(k.out$clust!=k.clust)]
      }
      j.length <- max(envi.cv)
      j.list[[j]] <- k.list[[which.max(lapply(X=c(1:j.length), FUN=function(X){k.list[[X]]$score}))]]$stack
      j.stack <- j.list[[j]]
      j.score <- max(unlist(lapply(X=c(1:j.length), FUN=function(X){k.list[[X]]$score})))
      j.vars <- rev(names(j.stack))
      j.v1 <- data.frame(t(data.frame(1:j.length))[0,])
      j.v2 <- data.frame(t(data.frame(j.vars)))
      j.v3 <- plyr::rbind.fill(j.v1, j.v2)
      colnames(j.v3) <- paste('var', 1:j.length, sep = '')
      j.v4 <- data.frame('algo'=n.algo, 'score'=j.score, j.v3)
      j.data <- rbind(j.data, j.v4)
    }
    i.t2 <- Sys.time(); i.time <- i.t2-i.t1
    t.list$time[which(t.list$algos==n.algo)] <- i.time
  }

  i.data <- plyr::ldply(.data=c(1:length(algos)), .fun=function(X){
    X.df <- data.frame('algo'=j.data$algo[X], 'score'=j.data$score[X],
                       'var'=data.frame(t(j.data[X, which(!colnames(j.data)%in%c('algo', 'score'))]))[,1])
    if(any(is.na(X.df$var))){X.df$var[which(is.na(X.df$var))] <- paste('NA', which(is.na(X.df$var)), sep='_')}
    return(X.df)})

  i.d2 <- plyr::ldply(.data=unique(i.data$var), .fun=function(X){
    data.frame('var'=X, 'score'=sum(i.data$score[which(i.data$var==X)]))
  })

  i.d2 <- i.d2[order(i.d2$score, decreasing=T),];  i.d2 <- i.d2[1:max(envi.cv),]
  i.stack <- c(envi2[[i.d2[1:max(envi.cv),'var'][!grepl('NA_',i.d2[1:max(envi.cv),'var'])]]])
  return(i.stack); options(digits=3)
}
