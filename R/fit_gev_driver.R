#===============================================================================
# Preliminary testing fitting GEV models.
# Note: this is antiquated by the "fit_many_gev_and_naveau_model.R" script.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#===============================================================================
# set up parameters
lnonstat <- vector('list',3); names(lnonstat) <- c('mu','sigma','xi')
lnonstat$mu    <- FALSE
lnonstat$sigma <- FALSE
lnonstat$xi    <- FALSE

# set parameter names, with an error check, that the nonstationary parameter
# combinations are supported
if(length(which(unlist(lnonstat)==TRUE))==0) {
  parnames <- c('mu','sigma','xi')
} else if(length(which(unlist(lnonstat)==TRUE))==1) {
  if(lnonstat$mu) {parnames <- c('mu0','mu1','sigma','xi')}
  else {print('ERROR - if only one nonstationary parameter, it must be location')}
} else if(length(which(unlist(lnonstat)==TRUE))==2) {
  if(lnonstat$mu & lnonstat$sigma) {parnames <- c('mu0','mu1','sigma0','sigma1','xi')}
  else {print('ERROR - if only two nonstationary parameters, they must be location and scale')}
} else if(length(which(unlist(lnonstat)==TRUE))==3) {
  parnames <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
} else {print('ERROR - unknown number of non-stationary parameters')}
#===============================================================================


#===============================================================================
# read temperature hindcast (1880-2016) and forecast (2016-2100) data
# yields time_forc, temperature_forc (temperature relative to 1901-2000, deg C)

library(ncdf4)
source('read_data_temperature.R')
#===============================================================================


#===============================================================================
# read and process the tide gauge data from Delfzijl

source('read_data_tidegauge.R')

data_calib <- lsl.max
auxiliary <- temperature_forc[1:length(data_calib)]
#===============================================================================


#===============================================================================
# set up some priors

# uniform priors
priors.unif <- vector('list', length(parnames)); names(priors.unif) <- parnames
for (p in 1:length(parnames)) {
  priors.unif[[p]] <- vector('list', 2); names(priors.unif[[p]]) <- c('lb','ub')
}
#===============================================================================


#===============================================================================
# fit MLE GEV forms

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)

source('likelihood_gev.R')

# get initial parameter estimates from MLE
gev.mle <- fevd(x=coredata(lsl.max))

##  Estimated parameters:
##     location        scale        shape
## 2849.54895053  442.37014103   -0.03305154
##  Standard Error Estimates:
##   location      scale      shape
## 43.20220492 31.56298254  0.07052645

# preliminary latin hypercube
bound.lower <- c(0, 0, -2)
bound.upper <- c(5000, 1000, 2)
niter.lhs <- 1e6
nparam.lhs <- length(bound.lower)

parameters.lhs0 <- randomLHS(n=niter.lhs, k=nparam.lhs)
parameters.lhs  <- parameters.lhs0
for (p in 1:nparam.lhs) {
  parameters.lhs[,p] <- parameters.lhs0[,p]*(bound.upper[p]-bound.lower[p]) + bound.lower[p]
}

llik.lhs <- rep(NA, niter.lhs)
lpost.lhs <- rep(NA, niter.lhs)
pb <- txtProgressBar(min=0,max=niter.lhs,initial=0,style=3)
for (i in 1:niter.lhs) {
  llik.lhs[i] <- log_like_gev(parameters=parameters.lhs[i,], parnames=parnames$gev3, data_calib=data_calib$lsl_max)
  lpost.lhs[i] <- log_post_gev(parameters=parameters.lhs[i,], parnames=parnames$gev3, data_calib=data_calib$lsl_max, model='gev3', priors=priors)
  setTxtProgressBar(pb, i)
}
close(pb)
lhs.sample <- data.frame(cbind(parameters.lhs, llik.lhs, lpost.lhs))
colnames(lhs.sample) <- c(parnames$gev3,'llike','lpost')

# sort lhs results from highest likelihood to lowest
lhs.llike <- lhs.sample[rev(order(llik.lhs)),]
lhs.lpost <- lhs.sample[rev(order(lpost.lhs)),]

# Differential evolution optimization
# (can only minimize, so use negative log-likelihood function)

# calculate BIC with each of the maximum likelihood sets (not that 'bestval' is
# from the negative log-likelihood, cancelling out the - in BIC)

NP.deoptim <- 100
niter.deoptim <- 100
F.deoptim <- 0.8
CR.deoptim <- 0.9

# all stationary
parnames <- c('mu','sigma','xi')
bound.lower <- c(0, 0, -3)
bound.upper <- c(8000, 4000, 3)
deoptim.gev3 <- DEoptim(neg_log_like_gev, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib)
names(deoptim.gev3$optim$bestmem) <- parnames
bic.gev3 <- 2*deoptim.gev3$optim$bestval + length(parnames)*log(length(data_calib))

# location nonstationary
parnames <- c('mu0','mu1','sigma','xi')
bound.lower <- c(0, -1000, 0, -3)
bound.upper <- c(8000, 1000, 4000, 3)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.gev4 <- DEoptim(neg_log_like_gev, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.gev4$optim$bestmem) <- parnames
bic.gev4 <- 2*deoptim.gev4$optim$bestval + length(parnames)*log(length(data_calib))

# location and scale nonstationary
parnames <- c('mu0','mu1','sigma0','sigma1','xi')
bound.lower <- c(0, -1000, 0, -200, -3)
bound.upper <- c(8000, 1000, 4000, 200, 3)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.gev5 <- DEoptim(neg_log_like_gev, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.gev5$optim$bestmem) <- parnames
bic.gev5 <- 2*deoptim.gev5$optim$bestval + length(parnames)*log(length(data_calib))

# location, scale and shape all nonstationary
parnames <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
bound.lower <- c(0, -1000, 0, -200, -3, -3)
bound.upper <- c(8000, 1000, 4000, 200, 3, 3)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.gev6 <- DEoptim(neg_log_like_gev, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.gev6$optim$bestmem) <- parnames
bic.gev6 <- 2*deoptim.gev6$optim$bestval + length(parnames)*log(length(data_calib))


#===============================================================================
# End
#===============================================================================
