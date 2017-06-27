#===============================================================================
# read tide gauge data for many locations to get a feel for what plausible
# values for the PP-GPD parameters may be
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

#rm(list=ls())
NP.deoptim000 <- 100      # number of DE population members (at least 10*[# parameters])
niter.deoptim000 <- 100   # number of DE iterations
#n_node000 <- 1            # number of CPUs to use
output.dir <- '../output/'
filename.processing <- '../output/preprocessing.RData'
#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/Users/axw322/codes/EVT/R')
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

print('reading and process data from Delfzijl and other European tide gauge stations.')
print(' (if you do not have this already done and [filename.processing] defined, this might take a while)...')

if(exists('filename.processing')) {load(filename.processing)
} else {source('processing_script.R')}

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
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
# for Dutch station
#===============================================================================
#

print('starting DE optimization for MLE PP-GPD parameters at Delfzijl, Netherlands...')

source('likelihood_ppgpd.R')

deoptim.delfzijl <- vector('list', nmodel); names(deoptim.delfzijl) <- types.of.model
for (i in 1:nmodel) {
  deoptim.delfzijl[[types.of.model[i]]] <- mat.or.vec(1, length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.delfzijl[[types.of.model[i]]]) <- 'delfzijl'
  colnames(deoptim.delfzijl[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.delfzijl <- mat.or.vec(1, nmodel)
if(nmodel > 1) {colnames(bic.delfzijl) <- types.of.model; rownames(bic.delfzijl) <- 'delfzijl'
} else {names(bic.delfzijl) <- 'delfzijl'}

# PP-GPD model fitting
for (gpd.type in types.of.gpd) {
#  if(gpd.type=='gpd3') {auxiliary <- NULL
#  } else {auxiliary <- trimmed_forcing(data_calib$gpd$year, time_forc, temperature_forc)$temperature}
# for now... (until add support for nonstationary)
  auxiliary <- NULL
  out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                       DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                       parnames=parnames_all[[gpd.type]], data_calib=data_calib, auxiliary=auxiliary)
  deoptim.delfzijl[[gpd.type]][1,] <- out.deoptim$optim$bestmem
  colnames(deoptim.delfzijl[[gpd.type]]) <- parnames_all[[gpd.type]]
  if(nmodel > 1) {
    bic.delfzijl[1, gpd.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_calib$gpd$counts_all))
  } else {
    bic.delfzijl <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_calib$gpd$counts_all))
  }
}

print('...done.')

# plot this survival function against the empirical one
if(FALSE) {
parameters <- deoptim.delfzijl$gpd3
x <- seq(0, 10000, 10)
nav.cdf <- naveau_cdf(x=x, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
plot(esf.levels, log10(esf.values), ylim=c(-2.5,0), xlim=c(0,6000), xlab='Level [cm]', ylab='Survival function [1-cdf]', yaxt='n')
axis(2, at=seq(-3, 0, by=1), label=parse(text=paste("10^", seq(-3,0), sep="")))
lines(x, log10(1-nav.cdf), col='red')
}

#
#===============================================================================
# fit MLE PP-GPD models for all tide gauge stations within +/-[some number of]
# degrees lat/lon of the Dutch station
#===============================================================================
#

print('starting DE optimization for MLE PP-GPD parameters for all European stations in set...')

deoptim.eur <- vector('list', nmodel); names(deoptim.eur) <- types.of.model
for (i in 1:nmodel) {
  deoptim.eur[[types.of.model[i]]] <- mat.or.vec(length(data_set), length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.eur[[types.of.model[i]]]) <- files.tg
  colnames(deoptim.eur[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.eur <- mat.or.vec(length(data_set), nmodel)
if(nmodel > 1) {colnames(bic.eur) <- types.of.model; rownames(bic.eur) <- files.tg
} else {names(bic.eur) <- files.tg}

for (dd in 1:length(data_set)) {
  print(paste('starting to calculate MLE PP-GPD parameters for European data set ',dd,' / ',length(data_set),sep=''))

  data_set[[dd]]$bic.deoptim <- rep(NA, nmodel)
  data_set[[dd]]$deoptim <- vector('list', nmodel)
  names(data_set[[dd]]$deoptim) <- types.of.model

  # PP-GPD model fitting
  for (gpd.type in types.of.gpd) {
    if(gpd.type=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gpd.type]], data_calib=data_europe[[dd]], auxiliary=auxiliary)
    deoptim.eur[[gpd.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[gpd.type]]) <- parnames_all[[gpd.type]]
    if(nmodel > 1) {
      bic.eur[dd, gpd.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_europe[[dd]]$gpd$counts_all))
    } else {
      bic.eur[dd] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gpd.type]])*log(length(data_set[[dd]]$gpd$counts_all))
    }
  }
}

# also include Delfzijl points in this fit/spread
mle.fits <- vector('list', nmodel); names(mle.fits) <- types.of.model
mle.fits$gpd3 <- rbind(deoptim.eur$gpd3, deoptim.delfzijl$gpd3)

print('...done.')

# plot a fit with the empirical survival function
if(FALSE){
dd <- 6
x.eur <- seq(from=0, to=6000, by=10)
parameters <- deoptim.eur$nav3[dd,]; cdf.nav <- naveau_cdf(x=x.eur, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
parameters <- deoptim.eur$gev3[dd,]; cdf.gev <- pevd(q=x.eur, loc=parameters[1], scale=parameters[2], shape=parameters[3])
esf.levels.eur <- data_set[[dd]]$lsl_max[order(data_set[[dd]]$lsl_max)]
esf.values.eur <- seq(from=length(data_set[[dd]]$lsl_max), to=1, by=-1)/(length(data_set[[dd]]$lsl_max)+1)
plot(esf.levels.eur, log10(esf.values.eur))
lines(x.eur, log10(1-cdf.nav), col='red')
lines(x.eur, log10(1-cdf.gev), col='blue')
}

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#


# TODO - HERE NOW
# TODO - HERE NOW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW!!
# TODO - HERE NOW


print('fitting prior distributions to the MLE parameters...')

# see where Delfzijl parameter fits lie with respect to the rest of them
# (make sure they are not some extreme outlier) (they are last row of each matrix)
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {
    hist(mle.fits[[model]][,p], xlab=parnames_all[[model]][p], main=model)
    lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
  }
}
}

# fit gamma and normal priors
# -> centered at the medians (or Delfzijl MLEs? no, medians, so it is general *prior* belief)
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if(exists('gamma.priors')) {rm(list=c('gamma.priors','normal.priors','uniform.priors'))}
gamma.priors <- c('lambda','lambda0','sigma','sigma0')
normal.priors <- c('xi','sigma1','xi0','xi1')
uniform.priors <- NULL

priors <- vector('list', nmodel); names(priors) <- types.of.model
for (model in types.of.model) {
  priors[[model]] <- vector('list', length(parnames_all[[model]])); names(priors[[model]]) <- parnames_all[[model]]
  for (par in parnames_all[[model]]) {
    priors[[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    if(!is.na(match(par, uniform.priors))) {
       names(priors[[model]][[par]]) <- c('type','shape','rate'); priors[[model]][[par]]$type <- 'uniform'
       priors[[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
       priors[[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
    } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
      names(priors[[model]][[par]]) <- c('type','shape','rate'); priors[[model]][[par]]$type <- 'gamma'
      priors[[model]][[par]]$rate <- median(mle.fits[[model]][,par]) / (0.5*(max(mle.fits[[model]][,par])-min(mle.fits[[model]][,par])))^2
      priors[[model]][[par]]$shape <- median(mle.fits[[model]][,par]) * priors[[model]][[par]]$rate
    } else if(!is.na(match(par, normal.priors))) {
      names(priors[[model]][[par]]) <- c('type','mean','sd'); priors[[model]][[par]]$type <- 'normal'
      priors[[model]][[par]]$mean <- median(mle.fits[[model]][,par])
      priors[[model]][[par]]$sd   <- 0.5*(max(mle.fits[[model]][,par])-min(mle.fits[[model]][,par]))
      #priors[[model]][[par]]$sd   <- sd(mle.fits[[model]][,par])
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
    range <- max(mle.fits[[model]][,p]) - min(mle.fits[[model]][,p])
    lower <- min(mle.fits[[model]][,p]) - 0.05*range
    upper <- max(mle.fits[[model]][,p]) + 0.05*range
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      range <- max(log(mle.fits[[model]][,p])) - min(log(mle.fits[[model]][,p]))
      lower <- min(log(mle.fits[[model]][,p])) - 0.05*range
      upper <- max(log(mle.fits[[model]][,p])) + 0.05*range
    }
    x.tmp <- seq(from=lower, to=upper, length.out=1000)
    if(priors[[model]][[p]]$type=='normal') {
      pdf.tmp <- dnorm(x=x.tmp, mean=priors[[model]][[p]]$mean, sd=priors[[model]][[p]]$sd)
    } else if(priors[[model]][[p]]$type=='gamma') {
      pdf.tmp <- dgamma(x=x.tmp, shape=priors[[model]][[p]]$shape, rate=priors[[model]][[p]]$rate)
    }
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      hist(log(mle.fits[[model]][,p]), xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(log(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p])), c(-1000,1000), type='l', col='red', lwd=2)
    } else {
      hist(mle.fits[[model]][,p], xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
    }
    lines(x.tmp, pdf.tmp, lwd=2, col='blue')
  }
}
}

#
#===============================================================================
# save priors and initial values (from Delfzijl DE optim.) and read later when
# calibrating with MCMC
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors <- paste(output.dir,'surge_priors_',appen,'_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_',appen,'_',today,'.rds', sep='')
filename.initvals <- paste(output.dir,'surge_initialvalues_',appen,'_',today,'.rds', sep='')
filename.everything <- paste(output.dir,'kitchen_sink_priors_',appen,'_',today,'.RData', sep='')

print(paste('saving priors and initial values as .rds files (',filename.priors,', ',filename.initvals,') to read and use later...',sep=''))

save.image(file=filename.everything)
saveRDS(priors, file=filename.priors)
saveRDS(mle.fits, file=filename.mles)
saveRDS(deoptim.delfzijl, file=filename.initvals)

print('...done.')

#
#===============================================================================
# End
#===============================================================================
#
