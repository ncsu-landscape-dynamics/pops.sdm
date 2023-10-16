#' @export

lc <- terra::rast('Q:\\Shared drives\\Data\\Raster\\USA\\landcover\\nlcd_2019_land_cover_l48_20210604_1s.tif')
vars <- c('built', 'decid', 'everg', 'trees', 'shrub', 'grass', 'pastr', 'cropl', 'culti', 'wetld')

lc.sum <- function(lc, var){
  if(var=='built'){vals <- c(21, 22, 23, 24)}
  if(var=='decid'){vals <- c(41, 43)}
  if(var=='everg'){vals <- c(42, 43)}
  if(var=='trees'){vals <- c(41, 42, 43)}
  if(var=='shrub'){vals <- c(52)}
  if(var=='grass'){vals <- c(71)}
  if(var=='pastr'){vals <- c(81)}
  if(var=='cropl'){vals <- c(82)}
  if(var=='culti'){vals <- c(81, 82)}
  if(var=='wetld'){vals <- c(90, 95)}

  lc.vals <- lc==vals

  if(var%in%c('built', 'decid', 'everg')){
    if(var=='built'){weights <- c(.25, .5, .75, 1)}
    if(var%in%c('decid', 'everg')){weights <- c(1, .5)}
    lc.vals <- (lc.vals)*weights
  }

  lc.var <- terra::app(lc.vals, fun='sum')

  terra::writeRaster(lc.var, paste('Q:\\Shared drives\\Data\\Raster\\USA\\landcover\\',
                                   paste('nlcd_2019_1s', var, paste(vals, collapse='_'), sep='_'),
                                   '.tif', sep=''))
}

i <- 1
for(i in i:length(vars)){lc.sum(lc, var=vars[[i]])}
