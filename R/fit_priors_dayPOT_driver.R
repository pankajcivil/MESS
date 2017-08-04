#===============================================================================
# fit_priors_dayPOT_driver.R
#
# Fit prior distributions for four candidate PP/GPD-based model structures, using
# all UHSLC data (research quality versions) at least 90 years, plus Sewell's
# Point (Norfolk) and Delfzijl. Makes for 30 stations total. (Balboa station is
# a part fo the UHSLC 'data_many' object)
#
# 1. get tide gauge data objects for all UHSLC database sites with > 90 years
# 2. get tide gauge data object for Delfzijl, the Netherlands
# 3. get ... for Norfolk, VA, USA (not in database)
# 4. get ... for Balboa, Panama (in the UHSLC network, so do not need to do anything here)
# 5. calculate maximum likelihood pp/gpd parameters, for each of four candidate
#    model structures, for each of the 30 sites
# 6. fit normal or gamma prior distributions to these 30 parameter sets, for each
#    model parameter within each of the four candidate model strcutures.
# 7. write this priors object to a file (rds) and save progress to revisit later
#    (rdata)
# 8. Also write a uniform priors object, using the same bounds as were imposed
#    on the MLE parameter search
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

#rm(list=ls())
NP.deoptim000 <- 100      # number of DE population members (at least 10*[# parameters])
niter.deoptim000 <- 100   # number of DE iterations
#n_node000 <- 1            # number of CPUs to use
output.dir <- '../output/'
filename.saveprogress <- '../output/fitting_priors.RData'
l.doprocessing <- 'FALSE'  # true if you need to run the processing
                          # false -> read in some previous RDS processing results,
                          # with the filenames defined below
filename.many <- '../data/tidegauge_processed_manystations_26Jul2017.rds'
filename.delfzijl <- '../data/tidegauge_processed_delfzijl_26Jul2017.rds'
filename.norfolk <- '../data/tidegauge_processed_norfolk_26Jul2017.rds'
# don't need filename.balboa because it is a part of 'data_many' object
#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/Users/tony/codes/EVT/R')
appen <- 'ppgpd'

#
#===============================================================================
# relevant libraries - do 'install.pacakges([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)
library(ncdf4)

#
#===============================================================================
# read and process data for temperature (auxiliary covariate for nonstationary
# parameters). yields: temperature_forc, time_forc
#===============================================================================
#

print('reading temperature data...')

source('read_data_temperature.R')

print('...done.')

#
#===============================================================================
# read and process data for tide gauge stations (or read previous processing)
#===============================================================================
#

print('reading/processing data from tide gauge stations...')
print(' (if you do not have this already done, this might take a while)...')

if(l.doprocessing) {
  source('processing_script.R')
} else {
  data_many <- readRDS(filename.many)
  data_delfzijl <- readRDS(filename.delfzijl)
  data_norfolk <- readRDS(filename.norfolk)
  data_balboa <- readRDS(filename.balboa)
}

# round them all up as one big data set
data_all <- data_many
data_all$delfzijl <- data_delfzijl
data_all$norfolk <- data_norfolk
# do not need to add Balboa because it is already part of 'data_many'

print('...done.')

#
#===============================================================================
# set up PP-GPD model parameters
#===============================================================================
#

print('setting up PP-GPD model parameters for DE optimization...')

source('parameter_setup_dayPOT.R')

print('...done.')

#
#===============================================================================
# parameters for DE optim (for maximum likelihood/minimum negative likelihoood)
#===============================================================================
#

#NP.deoptim <- 200
#niter.deoptim <- 200
NP.deoptim <- NP.deoptim000
niter.deoptim <- niter.deoptim000
F.deoptim <- 0.8
CR.deoptim <- 0.9

#
#===============================================================================
# fit MLE PP-GPD model parameters for each candidate model at each tide gauge
#===============================================================================
#

print('starting DE optimization for MLE PP-GPD parameters for all stations in set...')

# need the likelihood function
source('likelihood_ppgpd.R')

deoptim.all <- vector('list', nmodel); names(deoptim.all) <- types.of.model
for (i in 1:nmodel) {
  deoptim.all[[types.of.model[i]]] <- mat.or.vec(length(data_all), length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.all[[types.of.model[i]]]) <- names(data_all)
  colnames(deoptim.all[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.all <- mat.or.vec(length(data_all), nmodel)
if(nmodel > 1) {colnames(bic.all) <- types.of.model; rownames(bic.all) <- names(data_all)
} else {names(bic.all) <- names(data_all)}

for (dd in 1:length(data_all)) {
  print(paste('starting to calculate MLE PP-GPD parameters for tide gauge data set ',dd,' / ',length(data_all),sep=''))
  tbeg0 <- proc.time()
  data_all[[dd]]$bic.deoptim <- rep(NA, nmodel)
  data_all[[dd]]$deoptim <- vector('list', nmodel)
  names(data_all[[dd]]$deoptim) <- types.of.model
  # PP-GPD model fitting
  for (gpd.type in types.of.gpd) {
    print(paste('  - starting DE optimization for model ',gpd.type,'...', sep=''))
    tbeg <- proc.time()
    if(gpd.type=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_all[[dd]]$gpd$year, time_forc, temperature_forc)$temperature}
    # if tide gauge record starts before temperatures, clip it
    if(data_all[[dd]]$gpd$year[1] < time_forc[1]) {
      irem <- which(data_all[[dd]]$gpd$year < time_forc[1])
      auxiliary <- auxiliary[-irem]
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
      data_all[[dd]]$gpd$year <- data_all[[dd]]$gpd$year[-irem]
      data_all[[dd]]$gpd$counts <- data_all[[dd]]$gpd$counts[-irem]
      data_all[[dd]]$gpd$excesses <- data_all[[dd]]$gpd$excesses[-irem]
      data_all[[dd]]$gpd$time_length <- data_all[[dd]]$gpd$time_length[-irem]
      data_all[[dd]]$gpd$time_length_all <- sum(data_all[[dd]]$gpd$time_length)
      data_all[[dd]]$gpd$counts_all <- sum(unlist(data_all[[dd]]$gpd$counts), na.rm=TRUE)
      data_all[[dd]]$gpd$excesses_all <- unlist(data_all[[dd]]$gpd$excesses)[!is.na(unlist(data_all[[dd]]$gpd$excesses))]
    }
    out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gpd.type]], data_calib=data_all[[dd]]$gpd, auxiliary=auxiliary)
    deoptim.all[[gpd.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.all[[gpd.type]]) <- parnames_all[[gpd.type]]
    if(nmodel > 1) {
      bic.all[dd, gpd.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_all[[dd]]$gpd$counts_all))
    } else {
      bic.all[dd] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_set[[dd]]$gpd$counts_all))
    }
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  save.image(filename.saveprogress)
  tend0 <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))
}

print('...done.')

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#

print('fitting prior distributions to the MLE parameters...')

# fit gamma and normal priors
# -> centered at the medians
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if(exists('gamma.priors')) {rm(list=c('gamma.priors','normal.priors','uniform.priors'))}
gamma.priors <- c('lambda','lambda0','sigma','sigma0')
normal.priors <- c('lambda1','sigma1','xi','xi0','xi1')
uniform.priors <- NULL

priors_normalgamma <- vector('list', nmodel); names(priors_normalgamma) <- types.of.model
for (model in types.of.model) {
  priors_normalgamma[[model]] <- vector('list', length(parnames_all[[model]])); names(priors_normalgamma[[model]]) <- parnames_all[[model]]
  for (par in parnames_all[[model]]) {
    priors_normalgamma[[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    if(!is.na(match(par, uniform.priors))) {
       names(priors_normalgamma[[model]][[par]]) <- c('type','shape','rate'); priors_normalgamma[[model]][[par]]$type <- 'uniform'
       priors_normalgamma[[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
       priors_normalgamma[[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
    } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
      names(priors_normalgamma[[model]][[par]]) <- c('type','shape','rate'); priors_normalgamma[[model]][[par]]$type <- 'gamma'
      priors_normalgamma[[model]][[par]]$rate <- median(deoptim.all[[model]][,par]) / (0.5*(max(deoptim.all[[model]][,par])-min(deoptim.all[[model]][,par])))^2
      priors_normalgamma[[model]][[par]]$shape <- median(deoptim.all[[model]][,par]) * priors_normalgamma[[model]][[par]]$rate
    } else if(!is.na(match(par, normal.priors))) {
      names(priors_normalgamma[[model]][[par]]) <- c('type','mean','sd'); priors_normalgamma[[model]][[par]]$type <- 'normal'
      priors_normalgamma[[model]][[par]]$mean <- median(deoptim.all[[model]][,par])
      priors_normalgamma[[model]][[par]]$sd   <- 0.5*(max(deoptim.all[[model]][,par])-min(deoptim.all[[model]][,par]))
      #priors_normalgamma[[model]][[par]]$sd   <- sd(deoptim.all[[model]][,par])
    }
  }
}

print('...done.')

# plot priors and MLE histograms
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {
    range <- max(deoptim.all[[model]][,p]) - min(deoptim.all[[model]][,p])
    lower <- min(deoptim.all[[model]][,p]) - 0.05*range
    upper <- max(deoptim.all[[model]][,p]) + 0.05*range
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      range <- max(log(deoptim.all[[model]][,p])) - min(log(deoptim.all[[model]][,p]))
      lower <- min(log(deoptim.all[[model]][,p])) - 0.05*range
      upper <- max(log(deoptim.all[[model]][,p])) + 0.05*range
    }
    x.tmp <- seq(from=lower, to=upper, length.out=1000)
    if(priors_normalgamma[[model]][[p]]$type=='normal') {
      pdf.tmp <- dnorm(x=x.tmp, mean=priors_normalgamma[[model]][[p]]$mean, sd=priors_normalgamma[[model]][[p]]$sd)
    } else if(priors_normalgamma[[model]][[p]]$type=='gamma') {
      pdf.tmp <- dgamma(x=x.tmp, shape=priors_normalgamma[[model]][[p]]$shape, rate=priors_normalgamma[[model]][[p]]$rate)
    }
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      hist(log(deoptim.all[[model]][,p]), xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(log(c(deoptim.all[[model]][length(data_all),p],deoptim.all[[model]][length(data_all),p])), c(-1000,1000), type='l', col='red', lwd=2)
    } else {
      hist(deoptim.all[[model]][,p], xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(c(deoptim.all[[model]][length(data_all),p],deoptim.all[[model]][length(data_all),p]), c(-1000,1000), type='l', col='red', lwd=2)
    }
    lines(x.tmp, pdf.tmp, lwd=2, col='blue')
  }
}
}

#
#===============================================================================
# "fit" wide uniform priors (just using the bounds for the DE optim search)
#===============================================================================
#

print('fitting prior distributions to the uniform bounds for MLE search...')

# all parameters have uniform bounds, given by bound_lower_set and bound_upper_set
priors_uniform <- vector('list', nmodel); names(priors_uniform) <- types.of.model
for (model in types.of.model) {
  priors_uniform[[model]] <- vector('list', length(parnames_all[[model]])); names(priors_uniform[[model]]) <- parnames_all[[model]]
  for (par in parnames_all[[model]]) {
    priors_uniform[[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    names(priors_uniform[[model]][[par]]) <- c('type','lower','upper'); priors_uniform[[model]][[par]]$type <- 'uniform'
    priors_uniform[[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
    priors_uniform[[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
  }
}

print('...done.')

#
#===============================================================================
# save priors and initial values (from DE optim for Delfzijl, Balboa, and Norfolk)
# and read later when calibrating with MCMC
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors.normalgamma <- paste(output.dir,'surge_priors_normalgamma_',appen,'_',today,'.rds', sep='')
filename.priors.uniform <- paste(output.dir,'surge_priors_uniform_',appen,'_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_',appen,'_',today,'.rds', sep='')

print(paste('saving priors and DE optim output as .rds files to read and use later...',sep=''))

save.image(file=filename.saveprogress)
saveRDS(priors_normalgamma, file=filename.priors.normalgamma)
saveRDS(priors_uniform, file=filename.priors.uniform)
saveRDS(deoptim.all, file=filename.mles)


print('...done.')

#
#===============================================================================
# End
#===============================================================================
#
