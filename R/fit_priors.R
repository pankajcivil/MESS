#===============================================================================
# read tide gauge data for many locations to get a feel for what plausible
# values for the GEV and Naveau model (i) parameters may be
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

#rm(list=ls())
NP.deoptim000 <- 200      # number of DE population members (at least 10*[# parameters])
niter.deoptim000 <- 200   # number of DE iterations
#n_node000 <- 1            # number of CPUs to use
#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/Users/axw322/codes/EVT/R')
appen <- 'gev_nav'

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
# parameters)
# yields: temperature_forc, time_forc
#===============================================================================
#

print('reading temperature data...')

source('read_data_temperature.R')

print('...done.')

#
#===============================================================================
# read and process data for Dutch station (Delfzijl)
#===============================================================================
#

print('reading Delfzijl, Netherlands, tide gauge data...')

source('read_data_tidegauge_Delfzijl.R')

print('...done.')

#
#===============================================================================
# read and process data for set of European stations
#===============================================================================
#

print('reading lots of other European tide gauge data...')

filetype='csv'
dat.dir <- '~/codes/EVT/data/tide_gauge_Europe/'
files.tg <- list.files(path=dat.dir, pattern=filetype)

data_set <- vector('list', length(files.tg))
for (dd in 1:length(files.tg)) {
  names(data_set)[dd] <- substr(files.tg[dd], start=1, stop=7)
  data.tmp <- read.table(paste(dat.dir,files.tg[dd],sep=''), header=FALSE, sep=',')
  data_set[[dd]] <- vector('list', 5)
  names(data_set[[dd]]) <- c('year','month','day','hour','sl')
  data_set[[dd]]$year <- data.tmp$V1
  data_set[[dd]]$month <- data.tmp$V2
  data_set[[dd]]$day <- data.tmp$V3
  data_set[[dd]]$hour <- data.tmp$V4
  data_set[[dd]]$sl <- data.tmp$V5
}

# process data

# annual block maxima
for (dd in 1:length(data_set)) {
  data_set[[dd]]$sl_norm     <- rep(NA, length(data_set[[dd]]$year))
  data_set[[dd]]$year_unique <- unique(data_set[[dd]]$year)
  data_set[[dd]]$year_unique <- data_set[[dd]]$year_unique[order(data_set[[dd]]$year_unique)]
  data_set[[dd]]$lsl_max     <- rep(NA, length(data_set[[dd]]$year_unique))
  for (t in 1:length(data_set[[dd]]$year_unique)) {
    ind_this_year <- which(data_set[[dd]]$year==data_set[[dd]]$year_unique[t])
    data_set[[dd]]$sl_norm[ind_this_year] <- data_set[[dd]]$sl[ind_this_year] - mean(data_set[[dd]]$sl[ind_this_year])
    data_set[[dd]]$lsl_max[t] <- max(data_set[[dd]]$sl_norm[ind_this_year])
  }
}

# trim before 1850 (or whenever is time_forc[1]), which is when the temperature
# forcing starts. also make a note of how long each record is
nyear <- rep(NA, length(data_set))
for (dd in 1:length(data_set)) {
  if(data_set[[dd]]$year_unique[1] < time_forc[1]) {
    icut <- which(data_set[[dd]]$year_unique < time_forc[1])
    data_set[[dd]]$year_unique <- data_set[[dd]]$year_unique[-icut]
    data_set[[dd]]$lsl_max <- data_set[[dd]]$lsl_max[-icut]
  }
  nyear[dd] <- length(data_set[[dd]]$year_unique)
}

print('...done.')

#
#===============================================================================
# set up GEV and Naveau (i) model parameters
#===============================================================================
#

print('setting up GEV and Naveau-i model parameters for DE optimization...')

source('parameter_setup.R')

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

print('starting DE optimization for MLE GEV and Naveau-i parameters at Delfzijl, Netherlands...')

source('likelihood_gev.R')
source('likelihood_naveau.R')

deoptim.delfzijl <- vector('list', 8); names(deoptim.delfzijl) <- types.of.model
for (i in 1:length(types.of.model)) {
  deoptim.delfzijl[[types.of.model[i]]] <- mat.or.vec(1, length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.delfzijl[[types.of.model[i]]]) <- 'delfzijl'
  colnames(deoptim.delfzijl[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.delfzijl <- mat.or.vec(1, 8)
colnames(bic.delfzijl) <- types.of.model; rownames(bic.delfzijl) <- 'delfzijl'

# GEV model fitting
for (gev.type in types.of.gev) {
  if(gev.type=='gev3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  out.deoptim <- DEoptim(neg_log_like_gev, lower=bound_lower_set[[gev.type]], upper=bound_upper_set[[gev.type]],
                       DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                       parnames=parnames_all[[gev.type]], data_calib=data_calib$lsl_max, auxiliary=auxiliary)
  deoptim.delfzijl[[gev.type]][1,] <- out.deoptim$optim$bestmem
  colnames(deoptim.delfzijl[[gev.type]]) <- parnames_all[[gev.type]]
  bic.delfzijl[1, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gev.type]])*log(length(data_calib$lsl_max))
}

# Naveau (i) model fitting
for (nav.type in types.of.nav) {
  if(nav.type=='nav3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  out.deoptim <- DEoptim(neg_log_like_naveau, lower=bound_lower_set[[nav.type]], upper=bound_upper_set[[nav.type]],
                       DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                       parnames=parnames_all[[nav.type]], data_calib=data_calib$lsl_max, auxiliary=auxiliary)
  deoptim.delfzijl[[nav.type]][1,] <- out.deoptim$optim$bestmem
  colnames(deoptim.delfzijl[[nav.type]]) <- parnames_all[[nav.type]]
  bic.delfzijl[1, nav.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[nav.type]])*log(length(data_calib$lsl_max))
}

print('...done.')

# plot this survival function against the empirical one
if(FALSE) {
parameters <- deoptim.delfzijl$nav3
x <- seq(0, 10000, 10)
nav.cdf <- naveau_cdf(x=x, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
plot(esf.levels, log10(esf.values), ylim=c(-2.5,0), xlim=c(0,6000), xlab='Level [cm]', ylab='Survival function [1-cdf]', yaxt='n')
axis(2, at=seq(-3, 0, by=1), label=parse(text=paste("10^", seq(-3,0), sep="")))
lines(x, log10(1-nav.cdf), col='red')
}

#
#===============================================================================
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
# for all tide gauge stations within +/-15 degrees lat/lon of the Dutch station
#===============================================================================
#

print('starting DE optimization for MLE GEV and Naveau-i parameters for all European stations in set...')

deoptim.eur <- vector('list', 8); names(deoptim.eur) <- types.of.model
for (i in 1:length(types.of.model)) {
  deoptim.eur[[types.of.model[i]]] <- mat.or.vec(length(data_set), length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.eur[[types.of.model[i]]]) <- files.tg
  colnames(deoptim.eur[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.eur <- mat.or.vec(length(data_set), 8)
colnames(bic.eur) <- types.of.model; rownames(bic.eur) <- files.tg

for (dd in 1:length(data_set)) {

  data_set[[dd]]$bic.deoptim <- rep(NA, length(types.of.model))
  data_set[[dd]]$deoptim <- vector('list', length(types.of.model))
  names(data_set[[dd]]$deoptim) <- types.of.model

  # GEV model fitting
  for (gev.type in types.of.gev) {
    if(gev.type=='gev3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_gev, lower=bound_lower_set[[gev.type]], upper=bound_upper_set[[gev.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gev.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[gev.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[gev.type]]) <- parnames_all[[gev.type]]
    bic.eur[dd, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gev.type]])*log(length(data_set[[dd]]$lsl_max))
  }

  # Naveau (i) model fitting
  for (nav.type in types.of.nav) {
    if(nav.type=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_naveau, lower=bound_lower_set[[nav.type]], upper=bound_upper_set[[nav.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[nav.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[nav.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[nav.type]]) <- parnames_all[[nav.type]]
    bic.eur[dd, nav.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[nav.type]])*log(length(data_set[[dd]]$lsl_max))
  }
}

# also include Delfzijl points in this fit/spread
mle.fits <- vector('list', length(types.of.model)); names(mle.fits) <- types.of.model
mle.fits$gev3 <- rbind(deoptim.eur$gev3, deoptim.delfzijl$gev3)
mle.fits$gev4 <- rbind(deoptim.eur$gev4, deoptim.delfzijl$gev4)
mle.fits$gev5 <- rbind(deoptim.eur$gev5, deoptim.delfzijl$gev5)
mle.fits$gev6 <- rbind(deoptim.eur$gev6, deoptim.delfzijl$gev6)
mle.fits$nav3 <- rbind(deoptim.eur$nav3, deoptim.delfzijl$nav3)
mle.fits$nav4 <- rbind(deoptim.eur$nav4, deoptim.delfzijl$nav4)
mle.fits$nav5 <- rbind(deoptim.eur$nav5, deoptim.delfzijl$nav5)
mle.fits$nav6 <- rbind(deoptim.eur$nav6, deoptim.delfzijl$nav6)

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
rm(list=c('gamma.priors','normal.priors','uniform.priors'))
gamma.priors <- c('mu','mu0','kappa','kappa0','sigma','sigma0')
normal.priors <- c('xi','mu1','sigma1','xi0','xi1','kappa1')
uniform.priors <- NULL
# test fitting uniform priors
#uniform.priors <- c('mu','mu0','kappa','kappa0','sigma','sigma0','xi','mu1','sigma1','xi0','xi1','kappa1')

priors <- vector('list', length(types.of.model)); names(priors) <- types.of.model
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

# redo the kappas for Naveau, sampling from log(kappa) ('lkappa') instead
model <- 'nav3'; par <- 'kappa'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa'

model <- 'nav4'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

model <- 'nav5'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

model <- 'nav6'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

print('...done.')

# plot priors and MLE histograms
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {

# TODO - deal wtih lkappa and lkappa0 here

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
filename.priors <- paste('surge_priors_',appen,'_',today,'.rds', sep='')
filename.mles <- paste('surge_MLEs_',appen,'_',today,'.rds', sep='')
filename.initvals <- paste('surge_initialvalues_',appen,'_',today,'.rds', sep='')
filename.everything <- paste('kitchen_sink_priors_',appen,'_',today,'.RData', sep='')
filename.temperature <- paste('temperature_forcing_',today,'.csv', sep='')

print(paste('saving priors and initial values as .rds files (',filename.priors,', ',filename.initvals,') to read and use later...',sep=''))

save.image(file=filename.everything)
saveRDS(priors, file=filename.priors)
saveRDS(mle.fits, file=filename.mles)
saveRDS(deoptim.delfzijl, file=filename.initvals)
write.csv(x=cbind(time_forc, temperature_forc), file=filename.temperature, row.names=FALSE)

print('...done.')

#
#===============================================================================
# End
#===============================================================================
#
