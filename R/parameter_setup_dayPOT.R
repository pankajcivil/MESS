
types.of.gpd <- c('gpd3')
types.of.model <- c(types.of.gpd)
nmodel <- length(types.of.model)

# set up parameter names for each model
parnames_all <- vector('list', length(types.of.model))
parnames_all$gpd3 <- c('lambda','sigma','xi')

# set up parameter bounds for each model
bound_lower_set <- vector('list', length(types.of.model))
bound_lower_set$gpd3 <- c(0,0,-2)

bound_upper_set <- vector('list', length(types.of.model))
bound_upper_set$gpd3 <- c(1, 10000, 2)

#===============================================================================
# End
#===============================================================================
