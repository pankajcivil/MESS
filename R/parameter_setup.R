
types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

# if this is the first time we're doing this, set up so that we do the DE
# optimization from kappa itself...
if(!exists('priors')) {
  # sample log(kappa) or kappa?
  l_sample_log_kappa <- FALSE # FALSE -> do the transformatino later, for priors
# but if this is the actual calibration, we are sampling from log(kappa)
# Note that the parameter ranges (bound_lower/upper_set, below) don't matter for
# the MCMC calibration
} else if(exists('priors')) {
  l_sample_log_kappa <- TRUE
}

# set up parameter names for each model
parnames_all <- vector('list', length(types.of.model))
parnames_all$gev3 <- c('mu','sigma','xi')
parnames_all$gev4 <- c('mu0','mu1','sigma','xi')
parnames_all$gev5 <- c('mu0','mu1','sigma0','sigma1','xi')
parnames_all$gev6 <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
if(l_sample_log_kappa) {
  parnames_all$nav3 <- c('lkappa','sigma','xi')
  parnames_all$nav4 <- c('lkappa0','kappa1','sigma','xi')
  parnames_all$nav5 <- c('lkappa0','kappa1','sigma0','sigma1','xi')
  parnames_all$nav6 <- c('lkappa0','kappa1','sigma0','sigma1','xi0','xi1')
} else {
  parnames_all$nav3 <- c('kappa','sigma','xi')
  parnames_all$nav4 <- c('kappa0','kappa1','sigma','xi')
  parnames_all$nav5 <- c('kappa0','kappa1','sigma0','sigma1','xi')
  parnames_all$nav6 <- c('kappa0','kappa1','sigma0','sigma1','xi0','xi1')
}

# set up parameter bounds for each model
bound_lower_set <- vector('list', length(types.of.model))
bound_lower_set$gev3 <- c(0,0,-2)
bound_lower_set$gev4 <- c(0,-1000,0,-2)
bound_lower_set$gev5 <- c(0,-1000,0,-200,-2)
bound_lower_set$gev6 <- c(0,-1000,0,-200,-2,-3)
bound_lower_set$nav3 <- c(0,0,-2)
bound_lower_set$nav4 <- c(0,-100,0,-2)
bound_lower_set$nav5 <- c(0,-100, 0,-log(10),-2)
bound_lower_set$nav6 <- c(0,-100, 0,-log(10),-2,-3)
if(l_sample_log_kappa) {
  bound_lower_set$nav3 <- c(0,0,-2)                # in practice, log(kappa) > about 3
  bound_lower_set$nav4 <- c(0,-100,0,-2)
  bound_lower_set$nav5 <- c(0,-100, 0,-log(10),-2)
  bound_lower_set$nav6 <- c(0,-100, 0,-log(10),-2,-3)
}

bound_upper_set <- vector('list', length(types.of.model))
bound_upper_set$gev3 <- c(6000, 800, 2)
bound_upper_set$gev4 <- c(6000, 1000, 800, 2)
bound_upper_set$gev5 <- c(6000, 1000, 800, 200, 2)
bound_upper_set$gev6 <- c(6000, 1000, 800, 200, 2, 3)
bound_upper_set$nav3 <- c(1e5, 1000, 2)
bound_upper_set$nav4 <- c(1e5, 100, 1000, 2)
bound_upper_set$nav5 <- c(1e5, 100, log(1000), log(10), 2)
bound_upper_set$nav6 <- c(1e5, 100, log(1000), log(10), 2, 3)
if(l_sample_log_kappa) {
  bound_upper_set$nav3 <- c(log(1e5), 1000, 2)
  bound_upper_set$nav4 <- c(log(1e5), 100, 1000, 2)
  bound_upper_set$nav5 <- c(log(1e5), 100, log(1000), log(10), 2)
  bound_upper_set$nav6 <- c(log(1e5), 100, log(1000), log(10), 2, 3)
}


#===============================================================================
# End
#===============================================================================
