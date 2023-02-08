# require(biomod2)
# require(plyr)
# require(raster)
# require(terra)
#' @export

get_BestVars <- function(envi, pts){

  envi2 <- raster::stack(envi)
  pts.t <- which(pts$lyr1==1)
  pts.f <- which(pts$lyr1==0)
  p.max <- 1000000-length(pts.t)
  pts.r <- sample(x=pts.f, size=pmin(length(pts.f), p.max))
  pts.s <- pts[c(pts.t, pts.r)]
  myResp2 <- pts.s$lyr1
  myResp2[myResp2==0] <- NA # setting 'true absences' to NA
  myXY2 <- terra::crds(pts.s) # the XY coordinates of species data
  myName2 <- 'test'

  # algos <- list(c('SRE'), c('GLM'), c('RF'), c('MAXENT.Phillips'))
  # names(algos) <- c('SRE', 'GLM', 'RF', 'MAXENT.Phillips')
  algos <- list(c('SRE'), c('RF'), c('MAXENT'))
  names(algos) <- c('SRE', 'RF', 'MAXENT')

  i <- 1#; t.list <- list()
  k <- 1; k.list <- list()
  j <- 1; j.list <- list(); j.data <- data.frame()

  for(i in 1:length(algos)){
    i.algo <- algos[[i]]; n.algo <- names(algos)[i] #i.t1 <- Sys.time()
    for(j in 1:1){#reps of same algo
      for(k in 1:max(envi.cv$cluster)){ #reps through clusters
        if(k==1){k.names <- names(envi2); k.stack <-NULL}

        k.test <- plyr::ldply(.data=k.names,
                              .fun=function(X){
                                myExpl <- raster::stack(envi2[[X]], k.stack)
                                PA.df <- as.data.frame(myResp2); PA.df[is.na(PA.df)] <- FALSE
                                PA.fact <- sum(PA.df==F)/sum(myResp2, na.rm=T)
                                PA.df$myResp2[which(PA.df==F)[round((1:sum(myResp2, na.rm=T))*PA.fact)]] <- TRUE
                                PA.df$myResp2 <- as.logical(PA.df$myResp2)

                                myOptions <- biomod2::BIOMOD_ModelingOptions('GLM'=list(test='none'))
                                myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp2,
                                                                        expl.var = myExpl,
                                                                        resp.xy = myXY2,
                                                                        resp.name = myName2,
                                                                        PA.nb.rep = 1,
                                                                        PA.strategy = 'user.defined',
                                                                        PA.user.table = PA.df)

                                # 3. Computing the models
                                myModels <- biomod2::BIOMOD_Modeling(bm.format = myData,
                                                                     bm.options = myOptions,
                                                                     models = i.algo,
                                                                     nb.rep = 1, #number of runs
                                                                     data.split.perc = 100, #50,
                                                                     prevalence = NULL,
                                                                     var.import = 0,
                                                                     metric.eval = c('ROC', 'TSS'),
                                                                     save.output = T, # recommended to leave true
                                                                     scale.models = F, #experimental don't use
                                                                     do.full.models = F,
                                                                     modeling.id = paste(myName2,"Modeling",sep=""))

                                # Evaluation
                                myEval <- biomod2::get_evaluations(myModels) # get all models evaluation
                                eval.v <- myEval$calibration[myEval$metric.eval%in%c('ROC', 'TSS')]#eval.v <- c(myEval[c('ROC','TSS'),"Testing.data",,,])
                                return(mean(eval.v))
                              })
        row.names(k.test) <- k.names

        k.out <- data.frame(cbind(envi.cv$cluster[row.names(k.test)], k.test$V1))
        colnames(k.out) <- c('cluster', 'score')
        k.out <- k.out[order(k.out$cluster),]
        k.out[,3] <- as.integer(rank(-k.out$score))
        k.var <- row.names(k.out)[which.min(k.out$V3)]
        k.score <- k.out$score[which.max(k.out$score)]
        k.cluster <- k.out[k.var, 'cluster']
        k.envi <- envi2[[k.var]]
        k.stack <- raster::stack(k.envi, k.stack)
        k.list[[k]] <- list('score'=k.score, 'stack'=k.stack)

        # this reduces the number of variables in each cluster to the top 2 (redux var) most important
        if(k==1){redux <- 2; k.redux <- data.frame()
        for(i in 1:length(unique(k.out$cluster))){
          i.clust <- k.out[k.out$cluster==i,]
          i.redux <- i.clust[order(i.clust$score, decreasing=T),][1:redux,]
          k.redux <- rbind(k.redux, i.redux)
        }
        k.out <- k.redux
        }
        k.names <-  row.names(k.out)[which(k.out$cluster!=k.cluster)]
      }
      j.length <- max(envi.cv$cluster)
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
    # i.t2 <- Sys.time(); i.time <- i.t2-i.t1
    # t.list[[i]] <- i.time
  }

  i.data <- plyr::ldply(.data=c(1:length(algos)), .fun=function(X){
    X.df <- data.frame('algo'=j.data$algo[X], 'score'=j.data$score[X],
                       'var'=data.frame(t(j.data[X, which(!colnames(j.data)%in%c('algo', 'score'))]))[,1])
    if(any(is.na(X.df$var))){X.df$var[which(is.na(X.df$var))] <- paste('NA', which(is.na(X.df$var)), sep='_')}
    return(X.df)})

  i.d2 <- plyr::ldply(.data=unique(i.data$var), .fun=function(X){
    data.frame('var'=X, 'score'=sum(i.data$score[which(i.data$var==X)]))
  })

  i.d2 <- i.d2[order(i.d2$score, decreasing=T),];  i.d2 <- i.d2[1:max(envi.cv$cluster),]
  i.stack <- raster::stack(envi2[[i.d2[1:max(envi.cv$cluster),'var'][!grepl('NA_',i.d2[1:max(envi.cv$cluster),'var'])]]])
  return(i.stack)
}
