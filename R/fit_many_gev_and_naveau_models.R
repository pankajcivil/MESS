#===============================================================================
# read tide gauge data for many locations to get a feel for what plausible
# values for the GEV and Naveau model (i) parameters may be
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

# read data
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
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
#===============================================================================
#

NP.deoptim <- 200
niter.deoptim <- 200
F.deoptim <- 0.8
CR.deoptim <- 0.9

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

# set up parameter names for each model
parnames <- vector('list', length(types.of.model))
parnames$gev3 <- c('mu','sigma','xi')
parnames$gev4 <- c('mu0','mu1','sigma','xi')
parnames$gev5 <- c('mu0','mu1','sigma0','sigma1','xi')
parnames$gev6 <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
parnames$nav3 <- c('kappa','sigma','xi')
parnames$nav4 <- c('kappa0','kappa1','sigma','xi')
parnames$nav5 <- c('kappa0','kappa1','sigma0','sigma1','xi')
parnames$nav6 <- c('kappa0','kappa1','sigma0','sigma1','xi0','xi1')

# set up parameter bounds for each model
bound_lower_set <- vector('list', length(types.of.model))
bound_lower_set$gev3 <- c(0,0,-2)
bound_lower_set$gev4 <- c(0,-1000,0,-2)
bound_lower_set$gev5 <- c(0,-1000,0,-200,-2)
bound_lower_set$gev6 <- c(0,-1000,0,-200,-2,-3)
bound_lower_set$nav3 <- c(0,0,-10)
bound_lower_set$nav4 <- c(0,-100,0,-10)
bound_lower_set$nav5 <- c(0,-100,-100,-100,-10)
bound_lower_set$nav6 <- c(0,-100,-100,-100,-10,-10)

bound_upper_set <- vector('list', length(types.of.model))
bound_upper_set$gev3 <- c(6000, 800, 2)
bound_upper_set$gev4 <- c(6000, 1000, 800, 2)
bound_upper_set$gev5 <- c(6000, 1000, 800, 200, 2)
bound_upper_set$gev6 <- c(6000, 1000, 800, 200, 2, 3)
bound_upper_set$nav3 <- c(2e7, 1000, 10)
bound_upper_set$nav4 <- c(2e7, 100, 1000, 10)
bound_upper_set$nav5 <- c(2e7, 100, 100, 100, 10)
bound_upper_set$nav6 <- c(2e7, 100, 100, 100, 10, 10)

deoptim.eur <- vector('list', 8); names(deoptim.eur) <- types.of.model
for (i in 1:length(types.of.model)) {
  deoptim.eur[[types.of.model[i]]] <- mat.or.vec(length(data_set), length(parnames[[types.of.model[i]]]))
  rownames(deoptim.eur[[types.of.model[i]]]) <- files.tg
  colnames(deoptim.eur[[types.of.model[i]]]) <- parnames[[types.of.model[i]]]
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
                         parnames=parnames[[gev.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[gev.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[gev.type]]) <- parnames[[gev.type]]
    bic.eur[dd, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames[[gev.type]])*log(length(data_set[[dd]]$lsl_max))
  }

  # Naveau (i) model fitting
  for (nav.type in types.of.nav) {
    if(nav.type=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_naveau, lower=bound_lower_set[[nav.type]], upper=bound_upper_set[[nav.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames[[nav.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[nav.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[nav.type]]) <- parnames[[nav.type]]
    bic.eur[dd, nav.type] <- 2*out.deoptim$optim$bestval + length(parnames[[nav.type]])*log(length(data_set[[dd]]$lsl_max))
  }
}

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

# also include Delfzijl points in this fit/spread
mle.fits <- vector('list', length(types.of.model)); names(mle.fits) <- types.of.model
mle.fits$gev3 <- rbind(deoptim.eur$gev3, deoptim.gev3$optim$bestmem)
mle.fits$gev4 <- rbind(deoptim.eur$gev4, deoptim.gev4$optim$bestmem)
mle.fits$gev5 <- rbind(deoptim.eur$gev5, deoptim.gev5$optim$bestmem)
mle.fits$gev6 <- rbind(deoptim.eur$gev6, deoptim.gev6$optim$bestmem)
mle.fits$nav3 <- rbind(deoptim.eur$nav3, deoptim.nav3$optim$bestmem)
mle.fits$nav4 <- rbind(deoptim.eur$nav4, deoptim.nav4$optim$bestmem)
mle.fits$nav5 <- rbind(deoptim.eur$nav5, deoptim.nav5$optim$bestmem)
mle.fits$nav6 <- rbind(deoptim.eur$nav6, deoptim.nav6$optim$bestmem)

# see where Delfzijl parameter fits lie with respect to the rest of them
# (make sure they are not some extreme outlier) (they are last row of each matrix)
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames[[model]])) {
    hist(mle.fits[[model]][,p], xlab=parnames[[model]][p], main=model)
    lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
  }
}
}

# TODO - fit gamma and normal priors
# -> centered at the medians (or Delfzijl MLEs? no, medians, so it is general *prior* belief)
# -> with standard deviation equal to half the max-min range

# assign which parameters have which priors
gamma.priors <- c('mu','mu0','kappa','kappa0','sigma','sigma0')
normal.priors <- c('xi','mu1','sigma1','xi0','xi1','kappa1')

priors <- vector('list', length(types.of.model)); names(priors) <- types.of.model
for (model in types.of.model) {
  priors[[model]] <- vector('list', length(parnames[[model]])); names(priors[[model]]) <- parnames[[model]]
  for (par in parnames[[model]]) {
    priors[[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    if(!is.na(match(par, gamma.priors))) {
      names(priors[[model]][[par]]) <- c('type','shape','rate'); priors[[model]][[par]]$type <- 'gamma'
      # TODO - FIT THE PARAMETERS HERE
      # TODO
    } else if(~is.na(match(par, normal.priors))) {
      names(priors[[model]][[par]]) <- c('type','mean','sd'); priors[[model]][[par]]$type <- 'normal'
      # TODO - FIT THE PARAMETERS HERE
      # TODO
    }
  }
}

# are there any correlations among the parameters that shoudl be accounted for in the priors?
#e.g., cor(mle.fits[[model]])
# -> not really

#===============================================================================
# End
#===============================================================================
