#===============================================================================
# analysis_and_plotting.R
#
# Analysis and plotting for wong, klufas and keller pp/gpd analysis.
# This assumes you've downloaded the Github repo and are running from the
# 'R' directory.
#
# If any this is vague/confusing and you have questions, email me. Seriously.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#===============================================================================
# Some preliminaries

setwd('~/codes/EVT/R')

library(ncdf4)
library(extRemes)

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

# set useful directories -- assumes you are in the 'R' directory within the repo
# directory structure
plot.dir <- '../figures/'
output.dir <- '../output/'

# calibrated parameter sets (samples; all should be same size)
# these are the results from 'calibration_dayPOT-experiments_driver.R'
filename.norfolk.normalgamma <-  paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_normalgamma_27Jul2017.nc', sep='')
filename.norfolk.uniform <-      paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_uniform_27Jul2017.nc',     sep='')
filename.balboa.normalgamma <-   paste(output.dir,'calibratedParameters_ppgpd-experiments_balboa_normalgamma_28Jul2017.nc',  sep='')
filename.balboa.uniform <-       paste(output.dir,'calibratedParameters_ppgpd-experiments_balboa_uniform_28Jul2017.nc',      sep='')
filename.delfzijl.normalgamma <- paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_normalgamma_29Jul2017.nc',sep='')
filename.delfzijl.uniform <-     paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_uniform_29Jul2017.nc',    sep='')

filename.norfolk.data <- '../data/tidegauge_processed_norfolk_26Jul2017.rds'
filename.balboa.data <- '../data/tidegauge_processed_balboa_26Jul2017.rds'
filename.delfzijl.data <- '../data/tidegauge_processed_delfzijl_26Jul2017.rds'
#===============================================================================



#===============================================================================
#===============================================================================
# ANALYSIS
#===============================================================================
#===============================================================================




#===============================================================================
# read temperature data for projecting surge levels covarying with global mean
# surface temperature
#===============================================================================

# yields 'temperature_forc' and 'time_forc', between 1850 and 2100
# temperatures are relative to 1901-2000 mean temperature.

source('read_data_temperature.R')

# check the normalization
temperature.check <- mean(temperature_forc[which(time_forc==1901):which(time_forc==2000)])
print(paste('mean(temperature_forc between 1901 and 2000) is: ', temperature.check, sep=''))
print('(should be normalized to 0 during that period)')

#===============================================================================
# read parameter sets
#===============================================================================

# create list objects to store the parameters
site.names <- c('Balboa','Delfzijl','Norfolk'); nsites <- length(site.names)
data.lengths <- c(30,50,70,90,110,130); ndata.exp <- length(data.lengths)
data.experiment.names <- rep(NA, length(data.lengths)); for (dd in 1:length(data.lengths)) {data.experiment.names[dd] <- paste('y',data.lengths[dd],sep='')}
types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6'); nmodel <- length(types.of.gpd)

# define a generic list object that will have values for (i) each site,
# (ii) each data length experiment and (iii) each GPD model structure. (in that
# order)
list.init <- vector('list', nsites); names(list.init) <- site.names
for (site in site.names) {
  list.init[[site]] <- vector('list', length(data.lengths)); names(list.init[[site]]) <- data.experiment.names
  for (data.exp in data.experiment.names) {
    list.init[[site]][[data.exp]] <- vector('list', nmodel); names(list.init[[site]][[data.exp]]) <- types.of.gpd
  }
}

# initialize anything you would like to be defined for each site, for each data
# length and for each GPD model structure
gpd.parameters <- list.init
rl100 <- list.init

# index to store within the list level [[data.experiments]] which of the
# experiments corresponds to using all of the data. so you would use:
#    gpd.parameters$delfzijl[[all.data]]$gpd3
# for example to get at the calibrated gpd3 (stationary model) parameters for
# Delfzijl.
all.data <- rep(NA, nsites); names(all.data) <- site.names

# hard-coding here, even though it is bad practice. the different sites ahve
# different lengths of record, making things a bit complicated.

# need the thresholds for exceedances used
threshold <- vector('list', nsites); names(threshold) <- site.names

# parameter nmaes
parnames <- vector('list', nmodel); names(parnames) <- types.of.gpd

for (site in site.names) {
  if(site=='Norfolk') {
    ncdata <- nc_open(filename.norfolk.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd89.',model,sep='')))
      all.data[[site]] <- 4
      parnames[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$y30$gpd3)
    data.tmp <- readRDS(filename.norfolk.data); threshold[[site]] <- data.tmp$gpd$threshold
  } else if(site=='Balboa') {ncdata <- nc_open(filename.balboa.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd90.',model,sep='')))
      gpd.parameters[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd107.',model,sep='')))
      all.data[[site]] <- 5
      data.tmp <- readRDS(filename.balboa.data); threshold[[site]] <- data.tmp$gpd$threshold
    }
    if(nrow(gpd.parameters[[site]]$y30$gpd3) != n.ensemble) {print('ERROR - all sites ensembles must be the same size')}
  } else if(site=='Delfzijl') {ncdata <- nc_open(filename.delfzijl.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd90.',model,sep='')))
      gpd.parameters[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd110.',model,sep='')))
      gpd.parameters[[site]]$y130[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd137.',model,sep='')))
      all.data[[site]] <- 6
      data.tmp <- readRDS(filename.delfzijl.data); threshold[[site]] <- data.tmp$gpd$threshold
    }
    if(nrow(gpd.parameters[[site]]$y30$gpd3) != n.ensemble) {print('ERROR - all sites ensembles must be the same size')}
  } else {print('ERROR - unrecognized site name')}
}



#===============================================================================
# calculate return levels
#===============================================================================

TODO

# years to grab the reutrn levels
rl.years <- c(2000, 2016, 2065); nyears <- length(rl.years)
year.names <- rep(NA, length(rl.years)); for (year in 1:length(rl.years)) {year.names[year] <- paste('y',rl.years[year],sep='')}
temperature.years <- temperature_forc[which(time_forc == rl.years)]
print('you probably just got an error message - dont worry about it, probably')


# need the functions to project the PP/GPD parameters
source('likelihood_ppgpd.R')

# this is definitely coded like a neanderthal. can come back and code it more
# efficiently..?
for (site in site.names) {
  for (data.len in data.experiment.names[1:all.data[[site]]]) {
    for (model in types.of.gpd) {
      rl100[[site]][[data.len]][[model]] <- vector('list', nyears); names(rl100[[site]][[data.len]][[model]]) <- year.names
      for (year in year.names) {
        rl100[[site]][[data.len]][[model]][[year]] <- rep(NA, n.ensemble)
        print(paste('calculating return levels...',site,data.len,model,year,sep=' - '))
        tbeg <- proc.time()
        pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
        for (sow in 1:n.ensemble) {

          sigma <- gpd.parameters[[site]][[data.len]][[model]][sow,3]

          rl100[[site]][[data.len]][[model]][[year]][sow] <- rlevd(100, scale=)

# TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW

          rlevd(period, loc = 0, scale = 1, shape = 0, threshold = 0,
          type = c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull",
          "Exponential", "Beta", "Pareto"),
          npy = 365.25, rate = 0.01)



          setTxtProgressBar(pb, sow)
        }
        close(pb)
      }
    }
  }
}






#===============================================================================
# Bayesian model averaging ensemble
#===============================================================================

TODO - modify below as needed


bic <- vector('list', nexp); names(bic) <- gpd.experiments
llik.mod <- vector('list', nexp); names(llik.mod) <- gpd.experiments
llik.mod.all <- vector('list', nexp); names(llik.mod.all) <- gpd.experiments
lpost.mod.all <- vector('list', nexp); names(lpost.mod.all) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  bic[[gpd.exp]] <- rep(NA, nmodel); names(bic[[gpd.exp]]) <- types.of.model
  llik.mod[[gpd.exp]] <- rep(NA, nmodel); names(llik.mod[[gpd.exp]]) <- types.of.model
  llik.mod.all[[gpd.exp]] <- vector('list', nmodel); names(llik.mod.all[[gpd.exp]]) <- types.of.model
  lpost.mod.all[[gpd.exp]] <- vector('list', nmodel); names(lpost.mod.all[[gpd.exp]]) <- types.of.model
}

for (gpd.exp in gpd.experiments) {
  for (model in types.of.model) {
    if(substr(model,4,4)!='3') {auxiliary <- trimmed_forcing(data_calib[[gpd.exp]]$year, time_forc, temperature_forc)$temperature}
    lpri.tmp <- rep(NA, n.ensemble[[model]])
    llik.tmp <- rep(NA, n.ensemble[[model]])
    llik.mod.all[[gpd.exp]][[model]] <- rep(NA, n.ensemble[[model]])
    lpost.mod.all[[gpd.exp]][[model]] <- rep(NA, n.ensemble[[model]])
    for (i in 1:n.ensemble[[model]]) {
      lpri.tmp[i] <- log_prior_ppgpd(parameters=parameters[[gpd.exp]][[model]][i,], parnames=parnames_all[[model]], priors=priors.dayPOT, model=model)
      llik.tmp[i] <- log_like_ppgpd(parameters=parameters[[gpd.exp]][[model]][i,], parnames=parnames_all[[model]], data_calib=data_calib[[gpd.exp]], auxiliary=auxiliary)
      ndata <- sum(data_calib[[gpd.exp]]$counts)
      llik.mod.all[[gpd.exp]][[model]][i] <- llik.tmp[i]
      lpost.mod.all[[gpd.exp]][[model]][i] <- llik.tmp[i] + lpri.tmp[i]
    }
    imax <- which.max(llik.tmp)
    bic[[gpd.exp]][[model]] <- -2*max(llik.tmp) + length(parnames_all[[model]])*log(ndata)
    llik.mod[[gpd.exp]][[model]] <- mean(llik.tmp[is.finite(llik.tmp)])
  }
}

# a few different model averaging schemes
# note that wgt_lik is actual BMA
wgt_bic <- vector('list', nexp); names(wgt_bic) <- gpd.experiments
wgt_lik <- vector('list', nexp); names(wgt_lik) <- gpd.experiments
reg_wgt_bic <- vector('list', nexp); names(reg_wgt_bic) <- gpd.experiments
reg_wgt_lik <- vector('list', nexp); names(reg_wgt_lik) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  wgt_bic[[gpd.exp]] <- rep(NA, nmodel); names(wgt_bic[[gpd.exp]]) <- types.of.model
  wgt_lik[[gpd.exp]] <- rep(NA, nmodel); names(wgt_lik[[gpd.exp]]) <- types.of.model
  reg_wgt_bic[[gpd.exp]] <- rep(NA, nmodel); names(reg_wgt_bic[[gpd.exp]]) <- types.of.model
  reg_wgt_lik[[gpd.exp]] <- rep(NA, nmodel); names(reg_wgt_lik[[gpd.exp]]) <- types.of.model
}

llik.ref <- rep(NA, nexp); names(llik.ref) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  llik.ref[[gpd.exp]] <- median(llik.mod[[gpd.exp]])
  wgt_bic[[gpd.exp]] <- min(bic[[gpd.exp]])/bic[[gpd.exp]]; wgt_bic[[gpd.exp]] <- wgt_bic[[gpd.exp]]/sum(wgt_bic[[gpd.exp]])
  wgt_lik[[gpd.exp]] <- exp(llik.mod[[gpd.exp]] - llik.ref[[gpd.exp]])/sum(exp(llik.mod[[gpd.exp]] - llik.ref[[gpd.exp]]))
}
#===============================================================================




#===============================================================================
#===============================================================================
# PLOTTING
#===============================================================================
#===============================================================================





#
#===============================================================================
# SOM FIGURE S1 -- show the histograms of the MLE parameters and superimpose the
#                  prior distribution (normal or gamma)


#===============================================================================
#


#
#===============================================================================
# FIGURE 1 – Comparison of the empirical survival function calculated from the
#            observed tide gauge data at Delfzijl (red points) against the
#            modeled survival function in the ensemble median (black points) and
#            5-95% credible range (error bars) for models (a) ST: all parameters
#            stationary, (b) NS1: lambda non-stationary, (c) NS3: lambda and
#            sigma non-stationary, (d) NS3: all non-stationary and (e) the
#            BMA-weighted ensemble.


#===============================================================================
#


#
#===============================================================================
# FIGURE 2 - Projected distributions in 2065 of (a) the 20-year storm surge
#            return level, (b) the 50-year return level and (c) the 100-year
#            return level, relative to 2015, for each of the candidate models
#            and the BMA-weighted ensemble. The ST model would show no increase,
#            so we have plotted the actual surge level in 2015 for the ST model.


#===============================================================================
#


#
#===============================================================================
# FIGURE 3 - Distributions of calibrated return levels (rows) for varying time
#            lengths of data employed for each candidate model (columns).


#===============================================================================
#


#
#===============================================================================
# FIGURE 4 - Bayesian model averaging weights (equation (XX)) for the four
#            candidate models, using (a) 30 years of tide gauge data from
#            Delfzijl, (b) 50 years of data, (c) 70 years of data, (d) 90 years
#            of data, (e) 110 years of data and (f) 137 years of data. Higher
#            values imply better model-data match.


#===============================================================================
#




#===============================================================================
#===============================================================================
#===============================================================================
# Old plots below here
#===============================================================================
#===============================================================================
#===============================================================================







#===============================================================================
# preliminary plots of projections... what's going on?
#=====================================================

# stationary case
xtmp <- seq(from=0, to=10000, by=10)
mle.gev3 <- amcmc_out$gev3[[1]]$samples[which.max(amcmc_out$gev3[[1]]$log.p),]
mle.nav3 <- amcmc_out$nav3[[1]]$samples[which.max(amcmc_out$nav3[[1]]$log.p),]
pdf.nav3 <- naveau_pdf(x=xtmp, kappa=exp(mle.nav3[1]), sigma=mle.nav3[2], xi=mle.nav3[3])
pdf.gev3 <- devd(x=xtmp, loc=mle.gev3[1], scale=mle.gev3[2], shape=mle.gev3[3])

plot(xtmp, pdf.gev3, type='l', col='blue'); lines(xtmp, pdf.nav3, type='l', col='red')

# above should verify that the stationary case looks the same (check this below)

# gev/nav-4 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev4 <- amcmc_out$gev4[[1]]$samples[which.max(amcmc_out$gev4[[1]]$log.p),]
mle.nav4 <- amcmc_out$nav4[[1]]$samples[which.max(amcmc_out$nav4[[1]]$log.p),]
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
mle.gev5 <- amcmc_out$gev5[[1]]$samples[which.max(amcmc_out$gev5[[1]]$log.p),]
mle.nav5 <- amcmc_out$nav5[[1]]$samples[which.max(amcmc_out$nav5[[1]]$log.p),]
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
mle.gev6 <- amcmc_out$gev6[[1]]$samples[which.max(amcmc_out$gev6[[1]]$log.p),]
mle.nav6 <- amcmc_out$nav6[[1]]$samples[which.max(amcmc_out$nav6[[1]]$log.p),]
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
# FIGURE
#
# how does the 2000-year return level change throughout 2000-2100?
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

#
# The actual figure
#

pdf(paste(plot.dir,'2000yReturnPeriod_projections.pdf',sep=''),width=8,height=3.5,colormodel='cmyk')
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
#===============================================================================



#===============================================================================
# End
#===============================================================================
