#' @export
ptsmean <- function(runx){return(mean(runx$validation - actual_obs_runs$validation))}

#' @export
combine_points <- function(actual_obs, sampled_sd_locations, group, threshold){
com_obs <- rbind(actual_obs, sampled_sd_locations[sampled_sd_locations$sd_perc == threshold & sampled_sd_locations$group == group, ])
com_obs <- terra::vect(x=com_obs, geom=c("x", "y"))
return(com_obs)
}
