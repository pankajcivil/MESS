#===============================================================================
# read tide gauge data for many locations to get a feel for what plausible
# values for the GEV parameters may be
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

# fit MLE GEVs (gev3=all stationary (i.e., 3 free parameters), gev4=location
# parameter is nonstationary, gev5=location and scale nonstationary, gev6=
# all three nonstationary)

NP.deoptim <- 100
niter.deoptim <- 100
F.deoptim <- 0.8
CR.deoptim <- 0.9

types.of.gev <- c('gev3','gev4','gev5','gev6')

deoptim.eur <- vector('list', 4); names(deoptim.eur) <- types.of.gev
deoptim.eur$gev3 <- mat.or.vec(length(data_set), 3); rownames(deoptim.eur$gev3) <- files.tg
deoptim.eur$gev4 <- mat.or.vec(length(data_set), 4); rownames(deoptim.eur$gev4) <- files.tg
deoptim.eur$gev5 <- mat.or.vec(length(data_set), 5); rownames(deoptim.eur$gev5) <- files.tg
deoptim.eur$gev6 <- mat.or.vec(length(data_set), 6); rownames(deoptim.eur$gev6) <- files.tg
bic.eur <- mat.or.vec(length(data_set), 4)
colnames(bic.eur) <- types.of.gev; rownames(bic.eur) <- files.tg

for (dd in 1:length(data_set)) {

  data_set[[dd]]$bic.deoptim <- rep(NA, length(types.of.gev))
  data_set[[dd]]$deoptim <- vector('list', length(types.of.gev))
  names(data_set[[dd]]$deoptim) <- types.of.gev

  for (gev.type in types.of.gev) {

    if(gev.type=='gev3') {
      parnames <- c('mu','sigma','xi')
      bound_lower_set <- c(0,0,-2)
      bound_upper_set <- c(6000,800,2)
      auxiliary <- NULL
    } else if(gev.type=='gev4') {
      parnames <- c('mu0','mu1','sigma','xi')
      bound_lower_set <- c(0,-200,0,-2)
      bound_upper_set <- c(6000,200,800,2)
      auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature
    } else if(gev.type=='gev5') {
      parnames <- c('mu0','mu1','sigma0','sigma1','xi')
      bound_lower_set <- c(0,-200,0,-200,-2)
      bound_upper_set <- c(6000,200,800,200,2)
      auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature
    } else if(gev.type=='gev6') {
      parnames <- c('mu0','mu1','sigma0','sigma1','xi0','xi1')
      bound_lower_set <- c(0,-200,0,-200,-2,-2)
      bound_upper_set <- c(6000,200,800,200,2,2)
      auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature
    } else {print('ERROR - unrecognized gev.type')}

    out.deoptim <- DEoptim(neg_log_like_gev, lower=bound_lower_set, upper=bound_upper_set,
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames, data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[gev.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[gev.type]]) <- parnames
    bic.eur[dd, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames)*log(length(data_set[[dd]]$lsl_max))

  }

}


# TODO - fit KDE (or something?) to the many parameter estimates
# TODO - also include the Hook of Holland points in this fit/spread

#===============================================================================
# End
#===============================================================================
