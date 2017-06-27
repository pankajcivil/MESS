#===============================================================================
# parameter_setup_dayPOT.R
#
# set up parameters, parameter names, and ranges for daily peaks-over-threshold
# model
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6')
types.of.model <- c(types.of.gpd)
nmodel <- length(types.of.model)

# set up parameter names for each model
parnames_all <- vector('list', length(types.of.model))
parnames_all$gpd3 <- c('lambda','sigma','xi')
parnames_all$gpd4 <- c('lambda0','lambda1','sigma','xi')
parnames_all$gpd5 <- c('lambda0','lambda1','sigma0','sigma1','xi')
parnames_all$gpd6 <- c('lambda0','lambda1','sigma0','sigma1','xi0','xi1')

# set up parameter bounds for each model
bound_lower_set <- vector('list', length(types.of.model))
bound_lower_set$gpd3 <- c(0,0,-2)
bound_lower_set$gpd4 <- c(0,-1,0,-2)
bound_lower_set$gpd5 <- c(0,-1,0,-200,-2)
bound_lower_set$gpd6 <- c(0,-1,0,-200,-2,-3)

bound_upper_set <- vector('list', length(types.of.model))
bound_upper_set$gpd3 <- c(1,800,2)
bound_upper_set$gpd4 <- c(1,1,800,2)
bound_upper_set$gpd5 <- c(1,1,800,200,2)
bound_upper_set$gpd6 <- c(1,1,800,200,2,3)

#===============================================================================
# End
#===============================================================================
