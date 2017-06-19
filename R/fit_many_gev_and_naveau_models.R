#===============================================================================
# read tide gauge data for many locations to get a feel for what plausible
# values for the GEV and Naveau model (i) parameters may be
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

rm(list=ls())

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

source('read_data_temperature.R')

# function to trim temperature forcing to fit TG record unique years
# there might be missing years in TG record, so need to match each year and not
# just plop down an evenly spaced sequence
trimmed_forcing <- function(year_tidegauge, year_temperature, temperature) {
  output <- vector('list', 2); names(output) <- c('time','temperature')
  # check the beginning
  if(year_temperature[1] > year_tidegauge[1]) {print('ERROR - tide gauge record starts before temperature; add support for this situation')}
  # check the end
  if(max(year_temperature) < max(year_tidegauge)) {print('ERROR - tide gauge record ends after temperature; add support for this situation')}
  # match the indices of year_tidegauge within year_temperature
  imatch <- match(year_tidegauge, year_temperature)
  output$time <- year_temperature[imatch]
  output$temperature <- temperature[imatch]
  return(output)
}

#
#===============================================================================
# read and process data for Dutch station (Delfzijl)
#===============================================================================
#

source('read_data_tidegauge.R')
data_calib <- vector('list', 2); names(data_calib) <- c('year_unique','lsl_max')
data_calib$year_unique <- year.unique
data_calib$lsl_max <- lsl.max

#
#===============================================================================
# read and process data for set of European stations
#===============================================================================
#

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

#
#===============================================================================
# set up GEV and Naveau (i) model parameters
#===============================================================================
#

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

# set up parameter names for each model
parnames_all <- vector('list', length(types.of.model))
parnames_all$gev3 <- c('mu','sigma','xi')
parnames_all$gev4 <- c('mu0','mu1','sigma','xi')
parnames_all$gev5 <- c('mu0','mu1','sigma0','sigma1','xi')
parnames_all$gev6 <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
parnames_all$nav3 <- c('kappa','sigma','xi')
parnames_all$nav4 <- c('kappa0','kappa1','sigma','xi')
parnames_all$nav5 <- c('kappa0','kappa1','sigma0','sigma1','xi')
parnames_all$nav6 <- c('kappa0','kappa1','sigma0','sigma1','xi0','xi1')

# set up parameter bounds for each model
bound_lower_set <- vector('list', length(types.of.model))
bound_lower_set$gev3 <- c(0,0,-2)
bound_lower_set$gev4 <- c(0,-1000,0,-2)
bound_lower_set$gev5 <- c(0,-1000,0,-200,-2)
bound_lower_set$gev6 <- c(0,-1000,0,-200,-2,-3)
#bound_lower_set$nav3 <- c(0,0,-10)
#bound_lower_set$nav4 <- c(0,-100,0,-10)
#bound_lower_set$nav5 <- c(0,-100, 0,-100,-10)
#bound_lower_set$nav6 <- c(0,-100, 0,-100,-10,-10)
bound_lower_set$nav3 <- c(0,0,-2)
bound_lower_set$nav4 <- c(0,-100,0,-2)
bound_lower_set$nav5 <- c(0,-100, 0,-100,-2)
bound_lower_set$nav6 <- c(0,-100, 0,-100,-2,-3)

bound_upper_set <- vector('list', length(types.of.model))
bound_upper_set$gev3 <- c(6000, 800, 2)
bound_upper_set$gev4 <- c(6000, 1000, 800, 2)
bound_upper_set$gev5 <- c(6000, 1000, 800, 200, 2)
bound_upper_set$gev6 <- c(6000, 1000, 800, 200, 2, 3)
#bound_upper_set$nav3 <- c(2e7, 1000, 10)
#bound_upper_set$nav4 <- c(2e7, 100, 1000, 10)
#bound_upper_set$nav5 <- c(2e7, 100, 100, 100, 10)
#bound_upper_set$nav6 <- c(2e7, 100, 100, 100, 10, 10)
bound_upper_set$nav3 <- c(1e5, 1000, 2)
bound_upper_set$nav4 <- c(1e5, 100, 1000, 2)
bound_upper_set$nav5 <- c(1e5, 100, 100, 100, 2)
bound_upper_set$nav6 <- c(1e5, 100, 100, 100, 2, 3)

#
#===============================================================================
# parameters for DE optim (for maximum likelihood/minimum negative likelihoood)
#===============================================================================
#

NP.deoptim <- 200
niter.deoptim <- 200
F.deoptim <- 0.8
CR.deoptim <- 0.9

#
#===============================================================================
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
# for Dutch station
#
# Fun fact: This section is a bit of a catch-22: the prior distributions, which
# are defined in the 'likelihood_xxx.R' routines below, are determined by doing
# this initial MLE fitting. This was run once with the prior distributions left
# empty and only the likelihood functions defined, then the prior distibutions
# were fit as below, and built into the 'likelihood_xxx.R' routines.
#===============================================================================
#

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


### TODO TODO TODO - make sure you're using the right priors
### TODO TODO TODO - (uniform only for testing)


# plot priors and MLE histograms
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {
    range <- max(mle.fits[[model]][,p]) - min(mle.fits[[model]][,p])
    lower <- min(mle.fits[[model]][,p]) - 0.05*range
    upper <- max(mle.fits[[model]][,p]) + 0.05*range
    x.tmp <- seq(from=lower, to=upper, length.out=1000)
    if(priors[[model]][[p]]$type=='normal') {
      pdf.tmp <- dnorm(x=x.tmp, mean=priors[[model]][[p]]$mean, sd=priors[[model]][[p]]$sd)
    } else if(priors[[model]][[p]]$type=='gamma') {
      pdf.tmp <- dgamma(x=x.tmp, shape=priors[[model]][[p]]$shape, rate=priors[[model]][[p]]$rate)
    }
    hist(mle.fits[[model]][,p], xlab=parnames_all[[model]][p], main=model, freq=FALSE)
    lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
    lines(x.tmp, pdf.tmp, lwd=2, col='blue')
  }
}
}

# are there any correlations among the parameters that shoudl be accounted for in the priors?
#e.g., cor(mle.fits[[model]])
# -> not really

#
#===============================================================================
# set up and run PRELIMINARY Markov chain Monte Carlo (MCMC) calibration
#===============================================================================
#

# first, do a set of single-chain preliminary calibrations to get estimates of
# the jump covariance matrix

nnode_mcmc <- 1
gamma_mcmc <- 0.5
niter_mcmc <- 1e5
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_prelim <- vector('list', length(types.of.model)); names(amcmc_prelim) <- types.of.model

for (model in types.of.gev) {
  initial_parameters <- as.numeric(deoptim.delfzijl[[model]])
  if(model=='gev3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
  tend=proc.time()
}

for (model in types.of.nav) {
  initial_parameters <- as.numeric(deoptim.delfzijl[[model]])
  if('lkappa'  %in% parnames_all[[model]]) {initial_parameters[match('lkappa' ,parnames_all[[model]])] <- log(initial_parameters[match('lkappa' ,parnames_all[[model]])])}
  if('lkappa0' %in% parnames_all[[model]]) {initial_parameters[match('lkappa0',parnames_all[[model]])] <- log(initial_parameters[match('lkappa0',parnames_all[[model]])])}
  if(model=='nav3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
  tend=proc.time()
}

#
#===============================================================================
# set up and run PRODUCTION MCMC calibration
#===============================================================================
#

# then use these initial estimates of step_mcmc to launch the parallel chains
# (from amcmc_prelim[[model]]$cov.jump)

nnode_mcmc <- 4
gamma_mcmc <- 0.5
niter_mcmc <- 2e3
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_out <- vector('list', length(types.of.model)); names(amcmc_out) <- types.of.model

for (model in types.of.model) {

  print(paste('Starting calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  if (model %in% types.of.gev) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    if(model=='gev3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_gev, niter_mcmc, initial_parameters, n.chain=4, n.cpu=4,
				            scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
  } else if (model %in% types.of.nav) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    if('lkappa'  %in% parnames_all[[model]]) {initial_parameters[match('lkappa' ,parnames_all[[model]])] <- log(initial_parameters[match('lkappa' ,parnames_all[[model]])])}
    if('lkappa0' %in% parnames_all[[model]]) {initial_parameters[match('lkappa0',parnames_all[[model]])] <- log(initial_parameters[match('lkappa0',parnames_all[[model]])])}
    if(model=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_naveau, niter_mcmc, initial_parameters, n.chain=4, n.cpu=4,
				             scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                             gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                             parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                             priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
  } else {print('error - unknown model type')}
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}



# test plot
if (FALSE) {
model <- 'nav4'
par(mfrow=c(6,2))
for (p in 1:length(parnames_all[[model]])) {
    plot(amcmc_out[[model]]$samples[,p], type='l', ylab=parnames_all[[model]][p])
    hist(amcmc_out[[model]]$samples[round(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all[[model]][p], main='')
}
}


# Gelman and Rubin diagnostics - determine and chop off for burn-in
if(FALSE) {

# TODO

    # TODO - tidy up and figure out how to store, burn-in, etc.

    chain1=amcmc.par1[[1]]$samples
    chain2=amcmc.par1[[2]]$samples
    chain3=amcmc.par1[[3]]$samples
    chain4=amcmc.par1[[4]]$samples

  niter.test = seq(from=10000, to=nrow(chain1), by=5000)
  gr.test = rep(NA,length(niter.test))

  for (i in 1:length(niter.test)){
    mcmc1 = as.mcmc(chain1[1:niter.test[i],])
    mcmc2 = as.mcmc(chain2[1:niter.test[i],])
    mcmc3 = as.mcmc(chain3[1:niter.test[i],])
    mcmc4 = as.mcmc(chain4[1:niter.test[i],])
    mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2, mcmc3, mcmc4))
    gr.test[i] = gelman.diag(mcmc_chain_list)[2]
  }
}

# Heidelberger and Welch diagnostics

# TODO

# heidel.diag(as.mcmc(chain1), eps=0.1, pvalue=0.05)



#
#===============================================================================
# Once satisfied they are converged, write samples to netCDF file
#===============================================================================
#

parameters.posterior <- vector('list', length(types.of.model)); names(parameters.posterior) <- types.of.model
covjump.posterior <- vector('list', length(types.of.model)); names(covjump.posterior) <- types.of.model
for (model in types.of.model) {

 # TODO
 # TODO

  # for now, just using from amcmc_prelim[[model]]$samples
  parameters.posterior[[model]] <- amcmc_prelim[[model]]$samples[round(0.5*nrow(amcmc_prelim[[model]]$samples)):nrow(amcmc_prelim[[model]]$samples),]
  covjump.posterior[[model]] <- amcmc_prelim[[model]]$cov.jump

 # TODO
 # TODO

}

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape we expect, just transpose it
lmax=0
for (model in types.of.model) {for (i in 1:length(parnames_all[[model]])){lmax=max(lmax,nchar(parnames_all[[model]][i]))}}

## Name the file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
appen <- ''
filename.parameters <- paste('evt_models_calibratedParameters_',appen,'_',today,'.nc',sep='')

dim.parameters <- vector('list', length(types.of.model)); names(dim.parameters) <- types.of.model
dim.parnames   <- vector('list', length(types.of.model)); names(dim.parnames)   <- types.of.model
var.parameters <- vector('list', length(types.of.model)); names(var.parameters) <- types.of.model
var.parnames   <- vector('list', length(types.of.model)); names(var.parnames)   <- types.of.model
var.covjump    <- vector('list', length(types.of.model)); names(var.covjump)    <- types.of.model
dim.ensemble   <- vector('list', length(types.of.model)); names(dim.ensemble)   <- types.of.model
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.time <- ncdim_def('ntime', '', (time_forc), unlim=FALSE)
var.time <- ncvar_def('time', '', dim.time, -999)
var.temperature <- ncvar_def('temperature', 'Global mean temperature relative to 1901-2000', dim.time, -999)
for (model in types.of.model) {
  dim.parameters[[model]] <- ncdim_def(paste('n.parameters.',model,sep=''), '', 1:length(parnames_all[[model]]), unlim=FALSE)
  dim.ensemble[[model]]   <- ncdim_def(paste('n.ensemble.',model,sep=''), 'ensemble member', 1:nrow(parameters.posterior[[model]]), unlim=FALSE)
  var.parameters[[model]] <- ncvar_def(paste('parameters.',model,sep=''), '', list(dim.parameters[[model]],dim.ensemble[[model]]), -999)
  var.parnames[[model]]   <- ncvar_def(paste('parnames.',model,sep=''), '', list(dim.name,dim.parameters[[model]]), prec='char')
  var.covjump[[model]]    <- ncvar_def(paste('covjump.',model,sep=''), '', list(dim.parameters[[model]],dim.parameters[[model]]), -999)
}
output.to.file <- vector('list', length(var.parameters) + length(var.parnames) + length(var.covjump) + 2)
output.to.file[[1]] <- var.time
output.to.file[[2]] <- var.temperature
cnt <- 3
for (model in types.of.model) {
  output.to.file[[cnt]] <- var.parameters[[model]]
  output.to.file[[cnt+1]] <- var.parnames[[model]]
  output.to.file[[cnt+2]] <- var.covjump[[model]]
  cnt <- cnt + 3
}
#outnc <- nc_create(filename.parameters, list(var.parameters, var.parnames, var.covjump))
outnc <- nc_create(filename.parameters, output.to.file)
ncvar_put(outnc, var.time, time_forc)
ncvar_put(outnc, var.temperature, temperature_forc)
for (model in types.of.model) {
  ncvar_put(outnc, var.parameters[[model]], t(parameters.posterior[[model]]))
  ncvar_put(outnc, var.parnames[[model]], parnames_all[[model]])
  ncvar_put(outnc, var.covjump[[model]], covjump.posterior[[model]])
}
nc_close(outnc)

# link straight into flood risk
parameters <- parameters.posterior

# ... or...

#
#===============================================================================
# read a previous calibration output file, to use for flood risk
#===============================================================================
#

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

parameters <- vector('list', length(types.of.model)); names(parameters) <- types.of.model
parnames_all <- vector('list', length(types.of.model)); names(parnames_all) <- types.of.model
covjump <- vector('list', length(types.of.model)); names(covjump) <- types.of.model
ncdata <- nc_open('evt_models_calibratedParameters__13Jun2017.nc')
  time_forc <- ncvar_get(ncdata, 'time')
  temperature_forc <- ncvar_get(ncdata, 'temperature')
  for (model in types.of.model) {
    parameters[[model]]   <- t(ncvar_get(ncdata, paste('parameters.',model,sep='')))
    parnames_all[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    covjump[[model]]      <- ncvar_get(ncdata, paste('covjump.',model,sep=''))
  }
nc_close(ncdata)

#===============================================================================
# End
#===============================================================================




if(FALSE) {
#
#===============================================================================
# Bonus materials
#===============================================================================
#


#===============================================================================
# preliminary plots of projections... what's going on?
#=====================================================

# stationary case
xtmp <- seq(from=0, to=10000, by=10)
mle.gev3 <- amcmc_prelim$gev3$samples[which.max(amcmc_prelim$gev3$log.p),]
mle.nav3 <- amcmc_prelim$nav3$samples[which.max(amcmc_prelim$nav3$log.p),]
pdf.nav3 <- naveau_pdf(x=xtmp, kappa=exp(mle.nav3[1]), sigma=mle.nav3[2], xi=mle.nav3[3])
pdf.gev3 <- devd(x=xtmp, loc=mle.gev3[1], scale=mle.gev3[2], shape=mle.gev3[3])

plot(xtmp, pdf.gev3, type='l', col='blue'); lines(xtmp, pdf.nav3, type='l', col='red')

# above should verify that the stationary case looks the same (check this below)

# gev/nav-4 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev4 <- amcmc_prelim$gev4$samples[which.max(amcmc_prelim$gev4$log.p),]
mle.nav4 <- amcmc_prelim$nav4$samples[which.max(amcmc_prelim$nav4$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev4 <- project_gev(parameters=mle.gev4, parnames=parnames_all$gev4, auxiliary=temperature_proj)
par.nav4 <- project_naveau(parameters=mle.nav4, parnames=parnames_all$nav4, auxiliary=temperature_proj)
par.tmp <- par.gev4[which(time_proj==2020),]; pdf.gev4.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2040),]; pdf.gev4.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2060),]; pdf.gev4.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2080),]; pdf.gev4.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2100),]; pdf.gev4.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2020),]; pdf.nav4.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2040),]; pdf.nav4.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2060),]; pdf.nav4.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2080),]; pdf.nav4.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2100),]; pdf.nav4.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])


# gev/nav-5 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev5 <- amcmc_prelim$gev5$samples[which.max(amcmc_prelim$gev5$log.p),]
mle.nav5 <- amcmc_prelim$nav5$samples[which.max(amcmc_prelim$nav5$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev5 <- project_gev(parameters=mle.gev5, parnames=parnames_all$gev5, auxiliary=temperature_proj)
par.nav5 <- project_naveau(parameters=mle.nav5, parnames=parnames_all$nav5, auxiliary=temperature_proj)
par.tmp <- par.gev5[which(time_proj==2020),]; pdf.gev5.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2040),]; pdf.gev5.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2060),]; pdf.gev5.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2080),]; pdf.gev5.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2100),]; pdf.gev5.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2020),]; pdf.nav5.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2040),]; pdf.nav5.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2060),]; pdf.nav5.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2080),]; pdf.nav5.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2100),]; pdf.nav5.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])


# gev/nav-6 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev6 <- amcmc_prelim$gev6$samples[which.max(amcmc_prelim$gev6$log.p),]
mle.nav6 <- amcmc_prelim$nav6$samples[which.max(amcmc_prelim$nav6$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev6 <- project_gev(parameters=mle.gev6, parnames=parnames_all$gev6, auxiliary=temperature_proj)
par.nav6 <- project_naveau(parameters=mle.nav6, parnames=parnames_all$nav6, auxiliary=temperature_proj)
par.tmp <- par.gev6[which(time_proj==2020),]; pdf.gev6.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2040),]; pdf.gev6.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2060),]; pdf.gev6.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2080),]; pdf.gev6.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2100),]; pdf.gev6.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2020),]; pdf.nav6.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2040),]; pdf.nav6.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2060),]; pdf.nav6.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2080),]; pdf.nav6.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2100),]; pdf.nav6.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])

par(mfrow=c(3,2))
plot(xtmp/1000, pdf.gev4.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3), axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, location nonstationary')
  lines(xtmp/1000, pdf.gev4.2040, col='blue'); lines(xtmp/1000, pdf.gev4.2060, col='purple');
  lines(xtmp/1000, pdf.gev4.2080, col='red');  lines(xtmp/1000, pdf.gev4.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav4.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, lower tail nonstationary')
  lines(xtmp/1000, pdf.nav4.2040, col='blue'); lines(xtmp/1000, pdf.nav4.2060, col='purple');
  lines(xtmp/1000, pdf.nav4.2080, col='red');  lines(xtmp/1000, pdf.nav4.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.gev5.2020, type='l', xlim=c(0,10), ylim=c(0, 2e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, location and scale nonstationary')
  lines(xtmp/1000, pdf.gev5.2040, col='blue'); lines(xtmp/1000, pdf.gev5.2060, col='purple');
  lines(xtmp/1000, pdf.gev5.2080, col='red');  lines(xtmp/1000, pdf.gev5.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav5.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, lower tail and scale nonstationary')
  lines(xtmp/1000, pdf.nav5.2040, col='blue'); lines(xtmp/1000, pdf.nav5.2060, col='purple');
  lines(xtmp/1000, pdf.nav5.2080, col='red');  lines(xtmp/1000, pdf.nav5.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)

plot(xtmp/1000, pdf.gev6.2020, type='l', xlim=c(0,10), ylim=c(0, 3e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, all nonstationary')
  lines(xtmp/1000, pdf.gev6.2040, col='blue'); lines(xtmp/1000, pdf.gev6.2060, col='purple');
  lines(xtmp/1000, pdf.gev6.2080, col='red');  lines(xtmp/1000, pdf.gev6.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav6.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, all nonstationary')
  lines(xtmp/1000, pdf.nav6.2040, col='blue'); lines(xtmp/1000, pdf.nav6.2060, col='purple');
  lines(xtmp/1000, pdf.nav6.2080, col='red');  lines(xtmp/1000, pdf.nav6.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)

#===============================================================================



#===============================================================================
# a nice projections plot
#========================

# want the 1:2000 event, ensemble 5%, 50%, and 95% quantiles

time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]

protection.target <- 1/2000

# save the quantiles for the ensemble for each year
return_period_target <- vector('list', length(types.of.model)); names(return_period_target) <- types.of.model
for (model in types.of.model) {
  return_period_target[[model]] <- mat.or.vec(length(time_proj), 3)
  colnames(return_period_target[[model]]) <- c('q5','q50','q95')
}


for (model in types.of.gev) {
  print(paste('starting projections for ',model,' now...',sep=''))
  pb <- txtProgressBar(min=0,max=length(time_proj),initial=0,style=3)
  for (t in 1:length(time_proj)) {
    parameters_project <- t(sapply(1:nrow(parameters[[model]]), function(i) {project_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj[t])}))
    colnames(parameters_project) <- c('mu','sigma','xi')
    level_target <- as.numeric(sapply(1:nrow(parameters_project), function(i) {qevd(p=1-protection.target, loc=parameters_project[i,'mu'], scale=parameters_project[i,'sigma'], shape=parameters_project[i,'xi'])}))
    return_period_target[[model]][t,] <- quantile(level_target, c(.05, .5, .95))
    setTxtProgressBar(pb, t)
  }
  close(pb)
}
for (model in types.of.nav) {
  print(paste('starting projections for ',model,' now...',sep=''))
  pb <- txtProgressBar(min=0,max=length(time_proj),initial=0,style=3)
  for (t in 1:length(time_proj)) {
    parameters_project <- t(sapply(1:nrow(parameters[[model]]), function(i) {project_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj[t])}))
    colnames(parameters_project) <- c('kappa','sigma','xi')
    level_target <- as.numeric(sapply(1:nrow(parameters_project), function(i) {naveau_invcdf(q=1-protection.target, kappa=parameters_project[i,'kappa'], sigma=parameters_project[i,'sigma'], xi=parameters_project[i,'xi'])}))
    return_period_target[[model]][t,] <- quantile(level_target, c(.05, .5, .95), na.rm=TRUE)
    setTxtProgressBar(pb, t)
  }
  close(pb)
}

# convert from mm to m
for (model in types.of.model) {return_period_target[[model]] <- return_period_target[[model]]/1000}


# TODO
# TODO -- plot the projections as shaded areas (polygon)
#  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gmsl.ctrl.hind.95[itmp],rev(gmsl.ctrl.hind.05[itmp])),
#          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
# TODO

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

pdf('2000yReturnPeriod_projections.pdf',width=8,height=3.5,colormodel='cmyk')
par(mfrow=c(1,2))
par(mai=c(.65,.65,.20,.2))
model <- 'gev3'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4.5,14.5), xlab='', ylab='', xaxt='n', yaxt='n', col='black', xaxs='i')
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(.5,.5,.5,.5), border=NA)
mtext('1/2000 surge level [m]', side=2, line=2.3, cex=1);
mtext('Year', side=1, line=2, cex=1);
axis(1,seq(2000,2100,by=20), cex.axis=1)
axis(2,seq(5,13, by=1), label=c('5','','7','','9','','11','','13'), cex.axis=1)

model <- 'gev4'; ic <- 6
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

model <- 'gev5'; ic <- 9
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

legend(2000, 15, c('stationary','location nonstationary', 'location, scale nonstationary', '5-95% credible range'),
       lty=c(1,1,1,NA), pch=c(NA,NA,NA,15), lwd=c(2,2,2,10), cex=1, col=c('black', mycol.rgb[6], mycol.rgb[9], rgb(.5,.5,.5)), bty='n' )

# Too crazy ... need to check?
#model <- 'gev6'; ic <- 7
#lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
#      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
#polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
#        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

par(mai=c(.65,.65,.20,.2))
model <- 'nav3'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4.5,14.5), xlab='', ylab='', xaxt='n', yaxt='n', col='black', xaxs='i')

model <- 'nav6'; ic <- 14
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(.5,.5,.5,.5), border=NA)
mtext('1/2000 surge level [m]', side=2, line=2.3, cex=1);
mtext('Year', side=1, line=2, cex=1);
axis(1,seq(2000,2100,by=20), cex.axis=1)
axis(2,seq(5,13, by=1), label=c('5','','7','','9','','11','','13'), cex.axis=1)

model <- 'nav4'; ic <- 6
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

model <- 'nav5'; ic <- 9
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

legend(2000, 15, c('stationary','lower tail nonstationary', 'lower tail, scale nonstationary', 'both tails, scale nonstationary', '5-95% credible range'),
       lty=c(1,1,1,1,NA), pch=c(NA,NA,NA,NA,15), lwd=c(2,2,2,2,10), cex=1, col=c('black', mycol.rgb[6], mycol.rgb[9], mycol.rgb[14], rgb(.5,.5,.5)), bty='n' )

dev.off()



model <- 'gev4'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4,8), xlab='Year', ylab='1/2000 flood level [m]', main='GEV, location nonstationary')
lines(time_proj, return_period_target[[model]][,'q5'], type='l', lwd=2, lty=2)
lines(time_proj, return_period_target[[model]][,'q95'], type='l', lwd=2, lty=2)

model <- 'nav4'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4,8), xlab='Year', ylab='1/2000 flood level [m]', main='Naveau, location nonstationary')
lines(time_proj, return_period_target[[model]][,'q5'], type='l', lwd=2, lty=2)
lines(time_proj, return_period_target[[model]][,'q95'], type='l', lwd=2, lty=2)



#as.numeric(sapply(1:length(time_proj), function(t) {qevd(p=protection.target, loc=parameters_project[t,'mu'], scale=parameters_project[t,'sigma'], shape=parameters_project[t,'xi'])}))

#p_fail <- 1-sapply(1:nrow(ss.gev), function(t) {pgev(q=1000*H_effective, xi=ss.gev[t,'shape'], mu=ss.gev[t,'location'], beta=ss.gev[t,'scale']) })[,1]

#===============================================================================


#===============================================================================
# Thinning?
#==========

# not required, but you can if you feel like it (generally only storage limitations)

# use the preliminary chains to get an estimate of about how many iterations we
# need to thin the chains
acf_cutoff <- 0.05
niter_thin <- rep(0, length(types.of.model)); names(niter_thin) <- types.of.model
for (model in types.of.model) {
  for (p in 1:length(parnames_all[[model]])) {
    acf_tmp <- acf(x=amcmc_prelim[[model]]$samples[round(0.5*niter_mcmc):niter_mcmc,p], plot=FALSE, lag.max=round(0.4*niter_mcmc))
    niter_thin[[model]] <- max(niter_thin[[model]], acf_tmp$lag[which(acf_tmp$acf < acf_cutoff)[1]])
  }
}





# Have the 'posteriors' as log(likelihood) + log(prior); if we want the
# log-likelihood (e.g., for calculating BIC) need to subtract the log-prior out:
# log(likelihood) = log.p - log(prior)
model <- 'nav4'
log.prior.tmp <- log_prior_naveau(parameters=amcmc_out[[model]]$samples[which.max(amcmc_out[[model]]$log.p),],
                        parnames=parnames_all[[model]],
                        priors=priors,
                        model=model,
                        auxiliary=auxiliary)



#
#===============================================================================
#
}
